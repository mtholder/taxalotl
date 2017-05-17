#!/usr/bin/env python
# from __future__ import print_function

INP_TAXONOMY_DIRNAME = '__inputs__'
MISC_DIRNAME = '__misc__'

####################################################################################################
# Some data (to later be refactored
EUK_MICROBE_DIRNAME = '__other__'
TOP_PARTS = ('Archaea',
             'Bacteria',
             'Eukaryota/Metazoa',
             'Eukaryota/Fungi',
             'Eukaryota/plants',
             'Eukaryota/Archaeplastida',
             'Eukaryota/' + MISC_DIRNAME,
             'Eukaryota/' + EUK_MICROBE_DIRNAME,
             'Viruses',
             '__misc__',
             )

_part_list = [('Life', TOP_PARTS, ''),]
METAZOA_PARTS = ('Annelida',
                 'Arthropoda',
                 'Bryozoa',
                 'Chordata',
                 'Cnidaria',
                 'Ctenophora',
                 'Mollusca',
                 'Nematoda',
                 'Platyhelminthes',
                 'Porifera',
                 '__misc__',
                 )
METAZOA = 'Metazoa'
METAZOA_FRAG = 'Eukaryota/{}/'.format(METAZOA)
METAZOA_PARTS  = tuple([METAZOA_FRAG + i for i in METAZOA_PARTS])
_part_list.append((METAZOA, METAZOA_PARTS, METAZOA_FRAG))

ARTHROPODA = 'Arthropoda'
ARTHROPODA_PARTS = ('Malacostraca', 'Arachnida', 'Insecta')
ARTHROPODA_FRAG = METAZOA_FRAG + ARTHROPODA + '/'
ARTHROPODA_PARTS  = tuple([ARTHROPODA_FRAG + i for i in ARTHROPODA_PARTS])
_part_list.append((ARTHROPODA, ARTHROPODA_PARTS, ARTHROPODA_FRAG))

INSECTA = 'Insecta'
INSECTA_PARTS = ('Diptera', 'Coleoptera', 'Lepidoptera', 'Hymenoptera')
INSECTA_FRAG = ARTHROPODA_FRAG + INSECTA + '/'
INSECTA_PARTS  = tuple([INSECTA_FRAG + i for i in INSECTA_PARTS])
_part_list.append((INSECTA, INSECTA_PARTS, INSECTA_FRAG))

# Data above here, to be refactored at some point
####################################################################################################
# Code below
from peyotl import get_logger
_LOG = get_logger(__name__)
PARTS_BY_NAME = {}
PART_FRAG_BY_NAME = {}
for key, parts, frag in _part_list:
    PARTS_BY_NAME[key] = parts
    PART_FRAG_BY_NAME[key] = frag
PART_NAMES = list(PARTS_BY_NAME.keys())
PART_NAMES.sort()
PART_NAMES = tuple(PART_NAMES)


def _write_taxon(header, dict_to_write, id_order, dest_path)
    if not set_to_write:
        _LOG.info('No records need to be written to "{}"'.format(dest_path))
        return
    _LOG.info('Writing {} records to "{}"'.format(len(id_order), dest_path))
    pd = os.path.split(dest_path)[0]
    assure_dir_exists(pd)
    with codecs.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
        for i in id_order:
            outp.write(dict_to_write[i])


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

class TaxonFileOnlyPartitionElement(PartitionElement):
    def __init__(self, path_pref, fragment, path_suffix, roots):
        PartitionElement.__init__(self, path_pref, fragment, path_suffix, roots)
    def add_synonym(self, el_id, line):
        self.add(el_id, line)
    def write_synonyms(self, header):
        pass

class TaxAndSynFileOnlyPartitionElement(PartitionElement):
    def __init__(self, path_pref, fragment, path_suffix, roots):
        PartitionElement.__init__(self, path_pref, fragment, path_suffix, roots)
        fn = 'taxonomy.tsv'
        assert path_suffix.endswith(fn)
        self.syn_path_suffix = path_suffix[:-len(fn)] + 'synonyms.tsv'
        self.syn_path = os.path.join(path_pref, fragment, INP_TAXONOMY_DIRNAME, self.syn_path_suffix)
        self.syn_stored = {}
        self.syn_id_order = []

    def add_synonym(self, el_id, line):
        self.syn_stored[el_id] = line
        self.syn_id_order.append(el_id)
    def write_synonyms(self, header):
        _write_taxon(header, self.syn_stored, self.syn_id_order, self.syn_path)



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


def do_partition(parts_dir,
                 res_wrapper,
                 part_name,
                 part_keys,
                 par_frag,
                 pe_class,
                 taxon_filename,
                 master_map,
                 parse_and_partition_fn
                 ):
    _LOG.info('part_keys = {}'.format(part_keys))
    _LOG.info('par_frag = {}'.format(par_frag))
    path_suffix = os.path.join(res_wrapper.id, taxon_filename)
    remove_input = True
    if part_name == 'Life':
        remove_input = False
        inp_filepath = os.path.join(res_wrapper.unpacked_filepath, taxon_filename)
        misc_par = parts_dir
    else:
        misc_par = os.path.join(parts_dir, par_frag)
        inp_filepath = os.path.join(misc_par, INP_TAXONOMY_DIRNAME, path_suffix)


    mapping = [(k, master_map[k]) for k in part_keys if k in master_map]
    if not mapping:
        _LOG.info("No {} mapping for {}".format(res_wrapper.id, part_name))
        return
    partition_el = []
    for tag, roots in mapping:
        partition_el.append(pe_class(path_pref=parts_dir,
                                    fragment=tag,
                                    path_suffix=path_suffix,
                                    roots=roots))
    partition_el.append(pe_class(misc_par, MISC_DIRNAME, path_suffix, None))
    for part in partition_el:
        o = part.existing_output
        if o:
            m = 'Output for {} already exists at "{}"'
            raise RuntimeError(m.format(part.fragment, o))
    tup = parse_and_partition_fn(inp_filepath, partition_el)
    id_by_par, id_to_el, id_to_line, syn_by_id, roots_set, garbage_bin, header, syn_header = tup

    finish_partition_from_dict(id_by_par, id_to_el, id_to_line, garbage_bin)
    register_synonyms(syn_by_id, id_to_el, garbage_bin)
    for part in partition_el:
        part.write_lines(header)
        pr = [r for r in roots_set if id_to_el.get(r) is part]
        part.write_roots(pr)
    if remove_input:
        _LOG.info("leaving cruft at {}".format(inp_filepath))
        # os.unlink(inp_filepath)


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