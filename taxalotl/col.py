from __future__ import print_function

import codecs
import os
import time
from peyotl import (add_or_append_to_dict, get_logger, assure_dir_exists)
from taxalotl.interim_taxonomy_struct import InterimTaxonomyData

_LOG = get_logger(__name__)

TOP_PARTS = ('archaea', 'bacteria', 'misc', 'viruses',
             'eukaryotes/animals',
             'eukaryotes/fungi',
             'eukaryotes/plants',
             'eukaryotes/other',
             'eukaryotes/misc')

class PartitionElement(object):
    def __init__(self, path_pref, fragment, path_suffix, roots):
        self.path_pref = path_pref
        self.fragment = fragment
        self.path_suffix = path_suffix
        self.dest_path = os.path.join(path_pref, fragment, '__inputs__', path_suffix)
        self.roots = roots
        self.all_stored = {}
        self.id_order = []
    def add(self, el_id, line):
        self.all_stored[el_id] = line
        self.id_order.append(el_id)
    def write_lines(self, header):
        if not self.all_stored:
            _LOG.info('No records need to be written to "{}"'.format(self.dest_path))
            return
        _LOG.info('Writing {} records to "{}"'.format(len(self.id_order), self.dest_path))
        pd = os.path.split(self.dest_path)[0]
        assure_dir_exists(pd)
        with codecs.open(self.dest_path, 'w', encoding='utf-8') as outp:
            outp.write(header)
            for i in self.id_order:
                outp.write(self.all_stored[i])

def partition_col(parts_dir, wrapper):
    col_mapping = (('archaea', frozenset([33524792])),
                   ('bacteria', frozenset([33521420])),
                   ('viruses', frozenset([33521407])),
                   ('eukaryotes/animals', frozenset([33521288])),
                   ('eukaryotes/fungi', frozenset([33521351])),
                   ('eukaryotes/plants', frozenset([33521293])),
                   ('eukaryotes/other', frozenset([33521595, 33523363])))
    partition_el = []
    col_taxon_filename = 'taxa.txt'
    path_suffix = os.path.join(wrapper.id, col_taxon_filename)
    for tag, roots in col_mapping:
        partition_el.append(PartitionElement(path_pref=parts_dir,
                                             fragment=tag,
                                             path_suffix=path_suffix,
                                             roots=roots))
    partition_el.append(PartitionElement(parts_dir, 'misc', path_suffix, None))
    col_taxon_filepath = os.path.join(wrapper.unpacked_filepath, col_taxon_filename)
    _partition_col_by_root_id(col_taxon_filepath, partition_el)

def _partition_col_by_root_id(complete_fp, partition_el_list):
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
    id_to_line = {}
    id_by_par = {}
    syn_by_id = {}
    id_to_el = {}
    with codecs.open(complete_fp, 'rU', encoding='utf-8') as inp:
        iinp = iter(inp)
        header = iinp.next()
        prev_line = None
        vt = unicode('\x0b') # Do some lines have vertical tabs? Of course they do....
        istwo = unicode('\x1e')
        for n, line in enumerate(iinp):
            if not line.endswith('\n'):
                if prev_line:
                    prev_line = prev_line + line[:-1]
                else:
                    prev_line = line[:-1]
                continue
            elif prev_line:
                line = prev_line + line
                prev_line = ''
            ls = line.split('\t')
            if n % 1000 == 0:
                _LOG.info(' read line {}'.format(n))
            try:
                col_id, accept_id, par_id = ls[0], ls[4], ls[5]
                col_id = int(col_id)
                if accept_id:
                    accept_id = int(accept_id)
                    syn_by_id.setdefault(accept_id, []).append((col_id, line))
                else:
                    if col_id in roots_set:
                        match_l = [i[1] for i in by_roots if col_id in i[0]]
                        assert len(match_l) == 1
                        match_el = match_l[0]
                        id_to_el[col_id] = match_el
                        match_el.add(col_id, line)
                    else:
                        assert par_id
                        par_id = int(par_id)
                        match_el = id_to_el.get(par_id)
                        if match_el is not None:
                            id_to_el[col_id] = match_el
                            match_el.add(col_id, line)
                        else:
                            id_by_par.setdefault(par_id, []).append(col_id)
                            id_to_line[col_id] = line
            except:
                _LOG.exception("Exception parsing line {}:\n{}".format(n, line))
                raise
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
    for accept_id, i_l_list in syn_by_id.items():
        match_el = id_to_el.get(accept_id)
        if match_el is None:
            match_el = garbage_bin
        for col_id, line in i_l_list:
            match_el.add(col_id, line)
    for part in partition_el_list:
        part.write_lines(header)