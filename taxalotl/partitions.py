#!/usr/bin/env python
from __future__ import print_function

import copy
import os

from peyotl import get_logger, read_as_json

from taxalotl.tax_partition import (INP_TAXONOMY_DIRNAME,
                                    MISC_DIRNAME,
                                    GEN_MAPPING_FILENAME,
                                    get_taxon_partition,
                                    use_tax_partitions)

_LOG = get_logger(__name__)
_LIFE = 'Life'
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
BASE_PARTITIONS_DICT = {_LIFE: _x}
del _x
NAME_TO_PARTS_SUBSETS = {}
NAME_TO_PARENT_FRAGMENT = {}
NONTERMINAL_PART_NAMES = []
TERMINAL_PART_NAMES = []
PART_NAME_TO_FRAGMENT = {}


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
    a_list = list(res.alias_list)
    base_res = res.base_resource
    if base_res:
        a_list.extend(base_res.alias_list)
    poss_ids = [res.id] + a_list + [res.base_id]
    for k in poss_ids:
        if k in master_mapping:
            return master_mapping[k]
    m = 'No entry for ids {} found in "{}".'
    raise RuntimeError(m.format(', '.join(poss_ids), fp))


def _fill_parts_indices(d, par_frag):
    global NAME_TO_PARTS_SUBSETS, NAME_TO_PARENT_FRAGMENT, NONTERMINAL_PART_NAMES
    for k, subd in d.items():
        NAME_TO_PARTS_SUBSETS[k] = tuple(subd.keys())
        NAME_TO_PARENT_FRAGMENT[k] = par_frag
        if par_frag:
            cf = os.path.join(par_frag, k)
        else:
            cf = k
        PART_NAME_TO_FRAGMENT[k] = cf
        if subd:
            NONTERMINAL_PART_NAMES.append(k)
            _fill_parts_indices(subd, cf)
        elif k != MISC_DIRNAME:
            TERMINAL_PART_NAMES.append(k)


_fill_parts_indices(BASE_PARTITIONS_DICT, '')
PART_NAMES = list(NAME_TO_PARTS_SUBSETS.keys())
PART_NAMES.sort()
PART_NAMES = tuple(PART_NAMES)
PREORDER_PART_LIST = tuple(NONTERMINAL_PART_NAMES)
POSTORDER_PART_LIST = tuple(reversed(PREORDER_PART_LIST))
NONTERMINAL_PART_NAMES.sort()
NONTERMINAL_PART_NAMES = tuple(NONTERMINAL_PART_NAMES)
TERMINAL_PART_NAMES.sort()
TERMINAL_PART_NAMES = tuple(TERMINAL_PART_NAMES)


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
def iter_existing_tax_dirs(path_pref, res_id):
    suffix = os.path.join(INP_TAXONOMY_DIRNAME, res_id)
    misc_suffix = os.path.join(MISC_DIRNAME, INP_TAXONOMY_DIRNAME, res_id)
    for df in PART_NAME_TO_FRAGMENT.values():
        p = os.path.join(path_pref, df, suffix)
        if os.path.exists(p):
            yield p
        p = os.path.join(path_pref, df, misc_suffix)
        if os.path.exists(p):
            yield p


def has_any_partition_dirs(path_pref, res_id):
    assert path_pref
    for p in iter_existing_tax_dirs(path_pref, res_id):
        return True
    return False


def find_partition_dirs_for_taxonomy(path_pref, res_id):
    return [i for i in iter_existing_tax_dirs(path_pref, res_id)]


def get_fragment_from_part_name(parts_key):
    return PART_NAME_TO_FRAGMENT[parts_key]


def get_part_inp_taxdir(parts_dir, part_key, taxonomy_id):
    df = get_fragment_from_part_name(part_key)
    return os.path.join(parts_dir, df, INP_TAXONOMY_DIRNAME, taxonomy_id)


def get_par_and_par_misc_taxdir(parts_dir, part_key, taxonomy_id):
    df = get_fragment_from_part_name(part_key)
    par_df = os.path.split(df)[0]
    misc_df = os.path.join(par_df, MISC_DIRNAME)
    par_part_key = os.path.split(par_df)[0]
    pmtd = os.path.join(parts_dir, misc_df, INP_TAXONOMY_DIRNAME, taxonomy_id)
    return par_part_key, pmtd


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


def write_info_for_res(outstream, res, part_name_to_split):
    _LOG.debug('part_name_to_split = {}'.format(part_name_to_split))
    par_frag = NAME_TO_PARENT_FRAGMENT[part_name_to_split]
    _LOG.debug('par_frag = {}'.format(par_frag))
    if par_frag and not res.has_been_partitioned_for_fragment(par_frag):
        par_name = os.path.split(par_frag)[-1]
        outstream.write(
            '{} does not cover or has not been partitioned into {}\n'.format(res.id, par_name))
        return
    part_keys = NAME_TO_PARTS_SUBSETS[part_name_to_split]
    master_map = res.get_primary_partition_map()
    mapping = [(k, master_map[k]) for k in part_keys if k in master_map]
    if not mapping:
        outstream.write("No {} mapping for {}\n".format(res.id, part_name_to_split))
        return
    fragment = os.path.join(par_frag, part_name_to_split) if par_frag else part_name_to_split
    if not res.has_been_partitioned_for_fragment(fragment):
        outstream.write(
            '{} does not cover or has not been partitioned into {}\n'.format(res.id, fragment))
        return
    outstream.write('What can I say about {} at {} ? Great stuff...\n'.format(res.id, fragment))


def do_partition(res, part_name_to_split):
    """Partition a parent taxon into descendants and garbagebin (__misc__) dir

    :param res: a wrapper around the resource. Used for id, part_source_filepath, 
    :param part_name_to_split must be one of the hard-coded keys in NAME_TO_PARENT_FRAGMENT
    """
    _LOG.debug('part_name_to_split = {}'.format(part_name_to_split))
    par_frag = NAME_TO_PARENT_FRAGMENT[part_name_to_split]
    _LOG.debug('par_frag = {}'.format(par_frag))
    if par_frag and not res.has_been_partitioned_for_fragment(par_frag):
        par_name = os.path.split(par_frag)[-1]
        do_partition(res, par_name)
    part_keys = NAME_TO_PARTS_SUBSETS[part_name_to_split]
    master_map = res.get_primary_partition_map()
    mapping = [(k, master_map[k]) for k in part_keys if k in master_map]
    if not mapping:
        _LOG.info("No {} mapping for {}".format(res.id, part_name_to_split))
        return
    fragment = os.path.join(par_frag, part_name_to_split) if par_frag else part_name_to_split
    if res.has_been_partitioned_for_fragment(fragment):
        _LOG.info("Partition for fragment {} has already been done.".format(fragment))
        return
    tp = get_taxon_partition(res, fragment)
    if not par_frag:
        tp.external_input_fp = os.path.join(res.partition_source_dir, res.taxon_filename)
    tp.do_partition(mapping)


def check_partition(res, part_name_to_split):
    par_frag = NAME_TO_PARENT_FRAGMENT[part_name_to_split]
    part_keys = NAME_TO_PARTS_SUBSETS[part_name_to_split]
    master_map = res.get_primary_partition_map()
    fragment = os.path.join(par_frag, part_name_to_split) if par_frag else part_name_to_split
    if not res.has_been_partitioned_for_fragment(fragment):
        _LOG.info("Partition for fragment {} has not been done.".format(fragment))
        return True
    pop_subdirs = [k for k in part_keys if k in master_map]
    if not pop_subdirs:
        _LOG.info("No {} mapping for {}".format(res.id, part_name_to_split))
        return True
    with use_tax_partitions() as cache:
        misc = get_taxon_partition(res, fragment)
        cache.clear_without_flush(misc.cache_key)
        subs = [get_taxon_partition(res, os.path.join(fragment, k)) for k in pop_subdirs]
        unpart = get_taxon_partition(res, _LIFE)
        unpart.external_input_fp = os.path.join(res.partition_source_dir, res.taxon_filename)
        check_partition_union(fragment, misc, subs, unpart)


def check_partition_union(fragment, misc, subs, unpartitioned):
    slice_roots, slice_ids = misc._debug_validity_check()
    for p in subs:
        p_ids = p._debug_validity_check()[1]
        slice_ids.update(p_ids)
        _LOG.warn('{} IDs from {} bring total in {} up to {}'.format(len(p_ids), p.fragment,
                                                                     len(slice_ids), misc.fragment))
        for p_root_id, root_obj in p._roots.items():
            pr = root_obj['par_id']
            if pr not in slice_ids:
                slice_roots.add(p_root_id)
    unpartitioned._debug_check_subtree_ids(slice_roots, slice_ids)


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
