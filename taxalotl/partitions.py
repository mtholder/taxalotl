#!/usr/bin/env python
from __future__ import print_function

import copy
import os

from peyotl import get_logger, read_as_json

from taxalotl.tax_partition import (INP_TAXONOMY_DIRNAME,
                                    MISC_DIRNAME,
                                    GEN_MAPPING_FILENAME,
                                    get_taxon_partition)

_LOG = get_logger(__name__)

####################################################################################################
# Some data (to later be refactored
_x = {
    'Archaea': {},
    'Bacteria': {},
    'Eukaryota': {
        'Archaeplastida': {
            'Glaucophyta': {},
            'Rhodophyta': {},
            'Chloroplastida': {},
            MISC_DIRNAME: {},
        },
        'Fungi': {},
        'Haptophyta': {},
        'Metazoa': {
            'Annelida': {},
            'Arthropoda': {
                'Arachnida': {},
                'Malacostraca': {},
                'Insecta': {
                    'Diptera': {},
                    'Coleoptera': {},
                    'Lepidoptera': {},
                    'Hymenoptera': {},
                    MISC_DIRNAME: {},
                },
                MISC_DIRNAME: {},
            },
            'Bryozoa': {},
            'Chordata': {},
            'Cnidaria': {},
            'Ctenophora': {},
            'Mollusca': {},
            'Nematoda': {},
            'Platyhelminthes': {},
            'Porifera': {},
            MISC_DIRNAME: {},
        },
        'SAR': {},
        MISC_DIRNAME: {}
    },
    'Viruses': {},
    MISC_DIRNAME: {},
}
BASE_PARTITIONS_DICT = {"Life": _x}
del _x
PARTS_BY_NAME = {}
PART_FRAG_BY_NAME = {}
NONTERMINAL_PART_NAMES = []
PART_NAME_TO_DIRFRAG = {}


def get_inp_taxdir(parts_dir, frag, taxonomy_id):
    return os.path.join(parts_dir, frag, INP_TAXONOMY_DIRNAME, taxonomy_id)


def get_misc_inp_taxdir(parts_dir, frag, taxonomy_id):
    return os.path.join(parts_dir, frag, MISC_DIRNAME, INP_TAXONOMY_DIRNAME, taxonomy_id)


def get_all_taxdir_and_misc_uncles(parts_dir, frag, taxonomy_id):
    """Returns a list of dirs for this taxonomy_id starting at
    the `frag` directory, but also including the __misc__ subdirectories
     of is ancestral directories.
    This represents the set of directories that should hold the taxa
        for this fragment allowing for underclassification of taxa, but
        not misclassification into a non-ancestral group.
    """
    d = [get_inp_taxdir(parts_dir, frag, taxonomy_id)]
    if os.sep in frag:
        frag = os.path.split(frag)[0]
        while frag:
            d.append(get_inp_taxdir(parts_dir, os.path.join(frag, MISC_DIRNAME), taxonomy_id))
            if os.sep in frag:
                frag = os.path.split(frag)[0]
            else:
                break
    return d


def get_auto_gen_part_mapper(res):
    fp = os.path.join(res.partitioned_filepath, GEN_MAPPING_FILENAME)
    if not os.path.isfile(fp):
        m = 'Mapping file not found at "{}"\nRun the build-partitions-maps command.'
        raise RuntimeError(m.format(fp))
    master_mapping = read_as_json(fp)
    a_list = getattr(res, 'aliases', [])
    if a_list is None:
        a_list = []
    poss_ids = [res.id] + a_list + [res.base_id]
    for k in poss_ids:
        if k in master_mapping:
            return master_mapping[k]
    m = 'No entry for ids {} found in "{}".'
    raise RuntimeError(m.format(', '.join(poss_ids), fp))


def _fill_parts_indices(d, par_frag):
    global PARTS_BY_NAME, PART_FRAG_BY_NAME, NONTERMINAL_PART_NAMES
    for k, subd in d.items():
        PARTS_BY_NAME[k] = tuple(subd.keys())
        PART_FRAG_BY_NAME[k] = par_frag
        if par_frag:
            cf = os.path.join(par_frag, k)
        else:
            cf = k
        PART_NAME_TO_DIRFRAG[k] = cf
        if subd:
            NONTERMINAL_PART_NAMES.append(k)
            _fill_parts_indices(subd, cf)


_fill_parts_indices(BASE_PARTITIONS_DICT, '')
PART_NAMES = list(PARTS_BY_NAME.keys())
PART_NAMES.sort()
PART_NAMES = tuple(PART_NAMES)
PREORDER_PART_LIST = tuple(NONTERMINAL_PART_NAMES)
NONTERMINAL_PART_NAMES.sort()
NONTERMINAL_PART_NAMES = tuple(NONTERMINAL_PART_NAMES)


def _rec_populate(d_to_fill, key_to_filled_set):
    # _LOG.info('key_to_filled_set = {}'.format(key_to_filled_set))
    if MISC_DIRNAME in d_to_fill:
        del d_to_fill[MISC_DIRNAME]
    for key, subd in d_to_fill.items():
        filled_set = key_to_filled_set.get(key)
        if subd:
            _rec_populate(subd, key_to_filled_set)
            if not filled_set:
                cu = set()
                for v in subd.keys():
                    fsv = key_to_filled_set.get(v)
                    if fsv:
                        cu.update(fsv)
                if cu:
                    key_to_filled_set[key] = cu


def fill_empty_anc_of_mapping(mapping):
    # _LOG.info('mapping = {}'.format(mapping))
    s = copy.deepcopy(BASE_PARTITIONS_DICT)
    _rec_populate(s, mapping)


# Data above here, to be refactored at some point
####################################################################################################
# Code below


def find_partition_dirs_for_taxonomy(path_pref, res_id):
    suffix = os.sep + os.path.join('', INP_TAXONOMY_DIRNAME, res_id)
    return [i for i, sd, fl in os.walk(path_pref) if i.endswith(suffix)]


def get_relative_dir_for_partition(parts_key):
    return PART_NAME_TO_DIRFRAG[parts_key]


def get_part_inp_taxdir(parts_dir, part_key, taxonomy_id):
    df = get_relative_dir_for_partition(part_key)
    return os.path.join(parts_dir, df, INP_TAXONOMY_DIRNAME, taxonomy_id)


def get_par_and_par_misc_taxdir(parts_dir, part_key, taxonomy_id):
    df = get_relative_dir_for_partition(part_key)
    par_df = os.path.split(df)[0]
    misc_df = os.path.join(par_df, MISC_DIRNAME)
    par_part_key = os.path.split(par_df)[0]
    pmtd = os.path.join(parts_dir, misc_df, INP_TAXONOMY_DIRNAME, taxonomy_id)
    return par_part_key, pmtd


def get_root_ids_for_subset(tax_dir):
    rf = os.path.join(tax_dir, 'roots.txt')
    idset = set()
    if os.path.exists(rf):
        content = [int(i.strip()) for i in open(rf, 'r') if i.strip()]
        idset.update(content)
    return idset


def merge_and_write_taxon_partition_list(tp_list):
    if not tp_list:
        return
    fp_set = set()
    for tp in tp_list:
        fp = tp.taxon_fp
        if fp in fp_set:
            tp.append_write()
        else:
            tp.write()
            fp_set.add(fp)


def do_partition(res,
                 part_name,
                 part_keys,
                 par_frag,
                 master_map):
    """Partition a parent taxon into descendants and garbagebin (__misc__) dir

    :param res: a wrapper around the resource. Used for id, part_source_filepath, 
    :param part_name:
    :param part_keys:
    :param par_frag:
    :param master_map:
    :return:
    """
    fragment = os.path.join(par_frag, part_name)
    mapping = [(k, master_map[k]) for k in part_keys if k in master_map]
    if not mapping:
        _LOG.info("No {} mapping for {}".format(res.id, part_name))
        return
    tp = get_taxon_partition(res, fragment)
    tp.do_partition(mapping)
    tp.flush()


def get_inverse_misc_non_misc_dir_for_tax(inp_dir, tax_id):
    """ If given an unpartitioned dir, return (misc, False) otherwise (canonical, True)
    """
    misc_suffix = "/" + os.path.join(MISC_DIRNAME, INP_TAXONOMY_DIRNAME, tax_id)
    if inp_dir.endswith(misc_suffix):
        non_misc = inp_dir[:-len(misc_suffix)]
        return os.path.join(non_misc, INP_TAXONOMY_DIRNAME, tax_id), False
    non_misc_suffix = "/" + os.path.join(INP_TAXONOMY_DIRNAME, tax_id)
    assert inp_dir.endswith(non_misc_suffix)
    non_misc = inp_dir[:-len(non_misc_suffix)]
    return os.path.join(non_misc, MISC_DIRNAME, INP_TAXONOMY_DIRNAME, tax_id), True


'''
def partition_from_mapping(res, mapping, inp_dir, partition_parsing_fn, par_dir):
    """Returns a pair: a TaxonPartition element and the filepath to remove (or None)
    
    :param res:  resource wrapper
    :param mapping: list of pairs of taxon subdir name and root id set for this res
    :param inp_dir: input directory
    :param partition_parsing_fn: parsing function that returns TaxonPartition
    :param par_dir: output parent
    """
    misc_suffix = os.path.join(MISC_DIRNAME, INP_TAXONOMY_DIRNAME, res.id)
    in_misc = inp_dir.endswith(misc_suffix)
    taxon_filename = res.taxon_filename
    path_suffix = os.path.join(res.id, taxon_filename)
    inp_filepath = os.path.join(inp_dir, taxon_filename)
    top_life_dir = os.path.join(res.partitioned_filepath, 'Life')
    remove_input = (not in_misc) and (not top_life_dir == inp_dir)
    partition_el = []
    for tag, roots in mapping:
        pe = create_partition_element(path_pref=par_dir,
                                      fragment=tag,
                                      path_suffix=path_suffix,
                                      roots=roots,
                                      syn_filename=res.synonyms_filename)
        partition_el.append(pe)
    if in_misc:
        pe = create_partition_element(taxon_filepath=inp_filepath,
                                      roots=None,
                                      syn_filename=res.synonyms_filename)
    else:
        pe = create_partition_element(path_pref=par_dir,
                                      fragment=MISC_DIRNAME,
                                      path_suffix=path_suffix,
                                      roots=None,
                                      syn_filename=res.synonyms_filename)
    partition_el.append(pe)
    for part in partition_el:
        o = part.existing_output
        if o and not in_misc:
            m = 'Output for {} already exists at "{}"'
            raise RuntimeError(m.format(part.taxon_fp, o))
    if res.synonyms_filename:
        syn_file = os.path.join(os.path.split(inp_filepath)[0], res.synonyms_filename)
    else:
        syn_file = ''
    tp = TaxonPartition(sub=partition_el,
                        taxon_fp=inp_filepath,
                        syn_fp=syn_file,
                        parsing_func=partition_parsing_fn)
    to_remove_file = inp_filepath if remove_input and os.path.exists(inp_filepath) else None
    return tp, to_remove_file
'''
