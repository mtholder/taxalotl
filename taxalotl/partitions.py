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


def _append_taxon(dict_to_write, id_order, dest_path):
    if not dict_to_write:
        _LOG.info('No records need to be appended to "{}"'.format(dest_path))
        return
    _LOG.info('Appending {} records to "{}"'.format(len(id_order), dest_path))
    with codecs.open(dest_path, 'a', encoding='utf-8') as outp:
        if len(id_order) == len(dict_to_write):
            for i in id_order:
                outp.write(dict_to_write[i])
        else:
            for line in dict_to_write.values():
                outp.write(line)


def _write_taxon(header, dict_to_write, id_order, dest_path):
    _LOG.info('Writing header to "{}"'.format(dest_path))
    pd = os.path.split(dest_path)[0]
    assure_dir_exists(pd)
    with codecs.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
    _append_taxon(dict_to_write, id_order, dest_path)


def _write_taxon_list(header, record_list, dest_path):
    _LOG.info('Writing header to "{}"'.format(dest_path))
    pd = os.path.split(dest_path)[0]
    assure_dir_exists(pd)
    with codecs.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
    _append_taxon_list(record_list, dest_path)


def _append_taxon_list(record_list, dest_path):
    if not record_list:
        _LOG.info('No records need to be appended to "{}"'.format(dest_path))
        return
    _LOG.info('Appending {} records to "{}"'.format(len(record_list), dest_path))
    with codecs.open(dest_path, 'a', encoding='utf-8') as outp:
        for line in record_list:
            outp.write(line)


class PartitionElement(object):
    def __init__(self, path_pref, fragment, path_suffix, roots, dest_path):
        if dest_path:
            self.dest_path = dest_path
        else:
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

    def append_roots(self, root_ids):
        if not root_ids:
            _LOG.info('No root ids need to be apppended to "{}"'.format(self.roots_file))
            return
        sri = set(root_ids)
        with codecs.open(self.roots_file, 'r', encoding='utf-8') as outp:
            olri = [i.strip() for i in outp if i.strip()]
            oldset = set(olri)
            sri.update(oldset)
        self.write_roots(sri)

    def write_lines(self, header, syn_header=None):
        _write_taxon(header, self.all_stored, self.id_order, self.dest_path)
        self.write_synonyms(syn_header)

    def append_lines(self):
        _append_taxon(self.all_stored, self.id_order, self.dest_path)
        self.append_synonyms()

    @property
    def existing_output(self):
        if os.path.exists(self.dest_path):
            return self.dest_path
        return None

    def write_synonyms(self, header):
        pass

    def append_synonyms(self):
        pass

    def is_same_dest_pe(self, other):
        return self.dest_path == other.dest_path


class TaxonFileOnlyPartitionElement(PartitionElement):
    def __init__(self, path_pref, fragment, path_suffix, roots, dest_path):
        PartitionElement.__init__(self,
                                  path_pref,
                                  fragment,
                                  path_suffix,
                                  roots,
                                  dest_path=dest_path)

    def add_synonym(self, el_id, line):
        self.add(el_id, line)


class TaxAndSynFileOnlyPartitionElement(PartitionElement):
    def __init__(self, path_pref, fragment, path_suffix, roots, syn_filename, dest_path):
        PartitionElement.__init__(self,
                                  path_pref,
                                  fragment,
                                  path_suffix,
                                  roots,
                                  dest_path=dest_path)
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

    def append_synonyms(self):
        if self.id_less_syn:
            assert not self.syn_stored
            _append_taxon_list(self.id_less_syn, self.syn_path)
        else:
            _append_taxon(self.syn_stored, self.syn_id_order, self.syn_path)


def create_partition_element(path_pref=None,
                             fragment=None,
                             path_suffix=None,
                             roots=None,
                             syn_filename=None,
                             dest_path=None):
    if not syn_filename:
        return TaxonFileOnlyPartitionElement(path_pref,
                                             fragment,
                                             path_suffix,
                                             roots,
                                             dest_path=dest_path)
    return TaxAndSynFileOnlyPartitionElement(path_pref,
                                             fragment,
                                             path_suffix,
                                             roots,
                                             syn_filename,
                                             dest_path=dest_path)


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


def merge_and_write_taxon_partition_list(tp_list):
    if not tp_list:
        return
    dest_tp = tp_list[0]
    dest_tp.write()
    for another in tp_list[1:]:
        another.append_write()


class TaxonPartition(object):
    def __init__(self, taxon_fp, syn_fp, partition_el_list, parsing_func):
        self.partition_el = partition_el_list
        self.taxon_fp = taxon_fp
        self.syn_fp = syn_fp
        self.taxon_header = None
        self.syn_header = None
        t = separate_part_list(partition_el_list)
        self.roots_set, self.by_roots, self.garbage_bin = t
        self.id_to_line = {}
        self.id_by_par = {}
        self.syn_by_id = {}
        self.id_to_el = {}
        _LOG.info(repr(parsing_func))
        parsing_func(self)
        self.finish_partition_from_dict()
        self.register_synonyms()

    def finish_partition_from_dict(self):
        while True:
            par_id_matched = [p for p in self.id_by_par.keys() if p in self.id_to_el]
            if not par_id_matched:
                break
            for p in par_id_matched:
                id_list = self.id_by_par[p]
                match_el = self.id_to_el[p]
                for col_id in id_list:
                    self.id_to_el[col_id] = match_el
                    match_el.add(col_id, self.id_to_line[col_id])
                    del self.id_to_line[col_id]
                del self.id_by_par[p]
        if self.garbage_bin is not None:
            for col_id, line in self.id_to_line.items():
                self.garbage_bin.add(col_id, line)

    def register_synonyms(self):
        for accept_id, i_l_list in self.syn_by_id.items():
            match_el = self.id_to_el.get(accept_id)
            if match_el is None:
                match_el = self.garbage_bin
            for col_id, line in i_l_list:
                match_el.add_synonym(col_id, line)

    def write(self):
        for part in self.partition_el:
            part.write_lines(self.taxon_header, self.syn_header)
            pr = [r for r in self.roots_set if self.id_to_el.get(r) is part]
            part.write_roots(pr)

    def append_write(self):
        for part in self.partition_el:
            part.append_lines()
            pr = [r for r in self.roots_set if self.id_to_el.get(r) is part]
            part.append_write_roots(pr)


def do_new_separation(res,
                      new_par_dir,
                      inp_dir_list,
                      sub_dir_id_set_pairs_list):
    """Partition a parent taxon into descendants and garbagebin (__misc__) dir
    
    :param res: a wrapper around the resource. Used for id, part_source_filepath, 
    :param new_par_dir: parent of new dirs
    :param inp_dir_list: list of input dirs for this res type
    :param sub_dir_id_set_pairs_list: list of tuples the first element in each tuple
     is a name of a new subdir to be created. The second element in each tuple
     is the root IDs of this resource to be put into that dir.):
    """
    mapping = sub_dir_id_set_pairs_list
    if not mapping:
        _LOG.info("No {} mapping to separate {}".format(res.id, new_par_dir))
        return
    tp_list = []
    to_remove = []
    for inp_dir in inp_dir_list:
        tp, rm_file = _partition_from_mapping(res,
                                              mapping,
                                              inp_dir,
                                              partition_parsing_fn=res.partition_parsing_fn,
                                              par_dir=new_par_dir)
        tp_list.append(tp)
        to_remove.append(rm_file)
    merge_and_write_taxon_partition_list(tp_list)
    for to_remove_file in to_remove:
        if to_remove_file:
            _LOG.info("removing pre-partitioned file at {}".format(to_remove_file))
            os.unlink(to_remove_file)


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
    par_dir = os.path.join(res.partitioned_filepath, par_frag, part_name)
    taxon_filename = res.taxon_filename
    path_suffix = os.path.join(res.id, taxon_filename)
    if part_name == 'Life':
        inp_filepath = os.path.join(res.partition_source_filepath, taxon_filename)
    else:
        inp_filepath = os.path.join(par_dir, INP_TAXONOMY_DIRNAME, path_suffix)
    mapping = [(k, master_map[k]) for k in part_keys if k in master_map]
    if not mapping:
        _LOG.info("No {} mapping for {}".format(res.id, part_name))
        return
    tp, to_remove_file = _partition_from_mapping(res,
                                                 mapping,
                                                 os.path.split(inp_filepath)[0],
                                                 partition_parsing_fn=res.partition_parsing_fn,
                                                 par_dir=par_dir)
    tp.write()
    if to_remove_file is not None:
        _LOG.info("removing pre-partitioned file at {}".format(to_remove_file))
        os.unlink(to_remove_file)


def _partition_from_mapping(res, mapping, inp_dir, partition_parsing_fn, par_dir):
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
        pe = create_partition_element(dest_path=inp_dir,
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
            _LOG.info(m.format(part.fragment, o))
            return
    if res.synonyms_filename:
        syn_file = os.path.join(os.path.split(inp_filepath)[0], res.synonyms_filename)
    else:
        syn_file = ''
    tp = TaxonPartition(inp_filepath,
                        syn_file,
                        partition_el,
                        parsing_func=partition_parsing_fn)
    to_remove_file = inp_filepath if remove_input and os.path.exists(inp_filepath) else None
    return tp, to_remove_file
