#!/usr/bin/env python
from __future__ import print_function
from peyotl import get_logger, assure_dir_exists, read_as_json
import codecs
import copy
import os

_LOG = get_logger(__name__)
INP_TAXONOMY_DIRNAME = '__inputs__'
MISC_DIRNAME = '__misc__'
GEN_MAPPING_FILENAME = '__mapping__.json'

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
    #_LOG.info('key_to_filled_set = {}'.format(key_to_filled_set))
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
    #_LOG.info('mapping = {}'.format(mapping))
    s = copy.deepcopy(BASE_PARTITIONS_DICT)
    _rec_populate(s, mapping)


# Data above here, to be refactored at some point
####################################################################################################
# Code below


def _write_taxon(header, dict_to_write, id_order, dest_path):
    if not dict_to_write:
        _LOG.info('No records need to be written to "{}"'.format(dest_path))
        return
    _LOG.info('Writing {} records to "{}"'.format(len(id_order), dest_path))
    pd = os.path.split(dest_path)[0]
    assure_dir_exists(pd)
    with codecs.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
        if len(id_order) == len(dict_to_write):
            for i in id_order:
                outp.write(dict_to_write[i])
        else:
            for line in dict_to_write.values():
                outp.write(line)


def _write_taxon_list(header, record_list, dest_path):
    if not record_list:
        _LOG.info('No records need to be written to "{}"'.format(dest_path))
        return
    _LOG.info('Writing {} records to "{}"'.format(len(record_list), dest_path))
    pd = os.path.split(dest_path)[0]
    assure_dir_exists(pd)
    with codecs.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
        for line in record_list:
            outp.write(line)


class PartitionElement(object):
    def __init__(self, path_pref, fragment, path_suffix, roots):
        self.path_pref = path_pref
        self.fragment = fragment
        self.path_suffix = path_suffix
        self.dest_path = os.path.join(path_pref, fragment, INP_TAXONOMY_DIRNAME, path_suffix)
        self.roots = roots
        self.all_stored = {}
        self.id_order = []
        pd = os.path.split(self.dest_path)[0]
        self.roots_file = os.path.join(pd, 'roots.txt')

    def add(self, el_id, line):
        self.all_stored[el_id] = line
        self.id_order.append(el_id)

    def write_roots(self, root_ids):
        if not root_ids:
            _LOG.info('No root ids need to be written to "{}"'.format(self.roots_file))
            return
        _LOG.info('Writing {} root_ids to "{}"'.format(len(root_ids), self.roots_file))
        pd = os.path.split(self.roots_file)[0]
        assure_dir_exists(pd)
        with codecs.open(self.roots_file, 'w', encoding='utf-8') as outp:
            outp.write('\n'.join([str(i) for i in root_ids]))

    def write_lines(self, header, syn_header=None):
        _write_taxon(header, self.all_stored, self.id_order, self.dest_path)
        self.write_synonyms(syn_header)

    @property
    def existing_output(self):
        if os.path.exists(self.dest_path):
            return self.dest_path
        return None

    def write_synonyms(self, header):
        pass


class TaxonFileOnlyPartitionElement(PartitionElement):
    def __init__(self, path_pref, fragment, path_suffix, roots):
        PartitionElement.__init__(self, path_pref, fragment, path_suffix, roots)

    def add_synonym(self, el_id, line):
        self.add(el_id, line)


class TaxAndSynFileOnlyPartitionElement(PartitionElement):
    def __init__(self, path_pref, fragment, path_suffix, roots, syn_filename):
        PartitionElement.__init__(self, path_pref, fragment, path_suffix, roots)
        fn = 'taxonomy.tsv'
        assert path_suffix.endswith(fn)
        self.syn_path_suffix = path_suffix[:-len(fn)] + syn_filename
        self.syn_path = os.path.join(path_pref, fragment, INP_TAXONOMY_DIRNAME,
                                     self.syn_path_suffix)
        self.syn_stored = {}
        self.syn_id_order = []
        self.id_less_syn = []

    def add_synonym(self, el_id, line):
        if el_id:
            self.syn_stored[el_id] = line
            self.syn_id_order.append(el_id)
        else:
            self.id_less_syn.append(line)

    def write_synonyms(self, header):
        if self.id_less_syn:
            assert not self.syn_stored
            _write_taxon_list(header, self.id_less_syn, self.syn_path)
        else:
            _write_taxon(header, self.syn_stored, self.syn_id_order, self.syn_path)


def create_partition_element(path_pref, fragment, path_suffix, roots, syn_filename):
    if not syn_filename:
        return TaxonFileOnlyPartitionElement(path_pref, fragment, path_suffix, roots)
    return TaxAndSynFileOnlyPartitionElement(path_pref, fragment, path_suffix, roots, syn_filename)

def find_partition_dirs_for_taxonomy(path_pref, res_id):
    suffix = os.sep + os.path.join('', INP_TAXONOMY_DIRNAME, res_id)
    return [i for i, sd, fl in os.walk(path_pref) if i.endswith(suffix)]


def separate_part_list(partition_el_list):
    """Given a list of partition elements breaks the list baset on el.root for each el into:
        root_set = union of all el.roots
        by_roots list of (el.roots, el)
        garbage_bin the (max 1) element with el.roots == None
    """
    garbage_bin = None
    by_roots = []
    roots_set = set()
    for el in partition_el_list:
        if el.roots is None:
            assert garbage_bin is None
            garbage_bin = el
        else:
            by_roots.append((el.roots, el))
            roots_set.update(el.roots)
    return roots_set, by_roots, garbage_bin


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


def do_partition(res,
                 part_name,
                 part_keys,
                 par_frag,
                 master_map,
                 parse_and_partition_fn
                 ):
    """Partition a parent taxon into descendants and garbagebin (__misc__) dir
    
    :param res: a wrapper around the resource. Used for id, part_source_filepath, 
    :param part_name: 
    :param part_keys: 
    :param par_frag: 
    :param master_map: 
    :param parse_and_partition_fn: 
    :return: 
    """
    # _LOG.info('part_name = {}'.format(part_keys))
    # _LOG.info('part_keys = {}'.format(part_keys))
    # _LOG.info('par_frag = {}'.format(repr(par_frag)))
    par_dir = os.path.join(res.partitioned_filepath, par_frag)
    par_dir = os.path.join(par_dir, part_name)
    # _LOG.info('par_dir = {}'.format(repr(par_dir)))
    taxon_filename = res.taxon_filename
    # _LOG.info('taxon_filename = {}'.format(taxon_filename))
    path_suffix = os.path.join(res.id, taxon_filename)
    remove_input = True
    if part_name == 'Life':
        remove_input = False
        inp_filepath = os.path.join(res.partition_source_filepath, taxon_filename)
    else:
        inp_filepath = os.path.join(par_dir, INP_TAXONOMY_DIRNAME, path_suffix)

    mapping = [(k, master_map[k]) for k in part_keys if k in master_map]
    if not mapping:
        _LOG.info("No {} mapping for {}".format(res.id, part_name))
        return
    partition_el = []
    for tag, roots in mapping:
        pe = create_partition_element(path_pref=par_dir,
                                      fragment=tag,
                                      path_suffix=path_suffix,
                                      roots=roots,
                                      syn_filename=res.synonyms_filename)
        partition_el.append(pe)
    pe = create_partition_element(par_dir, MISC_DIRNAME, path_suffix, None, res.synonyms_filename)
    partition_el.append(pe)
    for part in partition_el:
        o = part.existing_output
        if o:
            m = 'Output for {} already exists at "{}"'
            _LOG.info(m.format(part.fragment, o))
            return
    if res.synonyms_filename:
        syn_file = os.path.join(os.path.split(inp_filepath)[0], res.synonyms_filename)
    else:
        syn_file = ''
    tup = parse_and_partition_fn(inp_filepath, syn_file, partition_el)
    id_by_par, id_to_el, id_to_line, syn_by_id, roots_set, garbage_bin, header, syn_header = tup

    finish_partition_from_dict(id_by_par, id_to_el, id_to_line, garbage_bin)
    register_synonyms(syn_by_id, id_to_el, garbage_bin)
    for part in partition_el:
        part.write_lines(header, syn_header)
        pr = [r for r in roots_set if id_to_el.get(r) is part]
        part.write_roots(pr)
    if remove_input and os.path.exists(inp_filepath):
        _LOG.info("removing unpartitioned taxon file at {}".format(inp_filepath))
        os.unlink(inp_filepath)


def register_synonyms(syn_by_id, id_to_el, garbage_bin):
    for accept_id, i_l_list in syn_by_id.items():
        match_el = id_to_el.get(accept_id)
        if match_el is None:
            match_el = garbage_bin
        for col_id, line in i_l_list:
            match_el.add_synonym(col_id, line)


def finish_partition_from_dict(id_by_par, id_to_el, id_to_line, garbage_bin):
    while True:
        par_id_matched = [p for p in id_by_par.keys() if p in id_to_el]
        if not par_id_matched:
            break
        for p in par_id_matched:
            id_list = id_by_par[p]
            match_el = id_to_el[p]
            for col_id in id_list:
                id_to_el[col_id] = match_el
                match_el.add(col_id, id_to_line[col_id])
                del id_to_line[col_id]
            del id_by_par[p]
    if garbage_bin is not None:
        for col_id, line in id_to_line.items():
            garbage_bin.add(col_id, line)
