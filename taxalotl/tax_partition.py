#!/usr/bin/env python
# from __future__ import print_function
import codecs
import os

from peyotl import get_logger, assure_dir_exists

_LOG = get_logger(__name__)


def _append_taxon(dict_to_write, id_order, dest_path):
    if not dict_to_write:
        _LOG.info('No records need to be appended to "{}"'.format(dest_path))
        return
    _LOG.info('Appending {} records to "{}"'.format(len(id_order), dest_path))
    with codecs.open(dest_path, 'a', encoding='utf-8') as outp:
        for i in id_order:
            el = dict_to_write.get(i)
            if el is not None:
                outp.write(el)
        if len(id_order) != len(dict_to_write):
            oset = frozenset(id_order)
            for key, line in dict_to_write.items():
                if key not in oset:
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


def _separate_part_list(partition_el_list):
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


class TaxonPartition(object):
    def __init__(self,
                 roots=None,
                 sub=None,
                 taxon_fp=None,
                 syn_fp=None,
                 parsing_func=None):
        self.id_order = []
        self.taxon_fp = taxon_fp
        self.syn_fp = syn_fp
        self.id_to_line = {}  # id -> line
        self.id_to_child_list = {}  # id -> list of child IDs
        self.syn_by_id = {}  # accepted_id -> list of synonym lines
        self.id_to_el = {}
        self.syn_id_order = []
        self.id_less_syn = []
        if sub is not None:
            assert roots is None
            self.partition_el = sub
            self.taxon_header = None
            self.syn_header = None
            t = _separate_part_list(self.partition_el)
            self.roots_set = t[0]  # union of all of the root sets of the partition elements
            self.by_roots = t[1]  # list of (ID set, PartitionElement) for each subset
            self.garbage_bin = t[3]  # the misc PartitionElement
            _LOG.info(repr(parsing_func))
            parsing_func(self)
            self._finish_partition_from_dict()
            self._register_synonyms()
            self.roots_file = None
        else:
            assert parsing_func is None
            assert roots is not None
            self.roots_set = roots
            if self.syn_fp is not None:
                fn = 'taxonomy.tsv'
                assert taxon_fp.endswith(fn)
            self.partition_el = None
            self.taxon_header = None
            self.syn_header = None
            self.by_roots = None
            self.garbage_bin = None
            pd = os.path.split(self.taxon_fp)[0]
            self.roots_file = os.path.join(pd, 'roots.txt')

    def add_taxon(self, uid, par_id, line):
        if uid not in self.id_to_line:
            self.id_to_line[uid] = line
            self.id_order.append(uid)
        self.id_to_child_list.setdefault(par_id, []).append(uid)

    def write_roots(self, root_ids):
        if self.garbage_bin is not None:
            self.garbage_bin.write_roots(root_ids)
            return
        if not root_ids:
            _LOG.info('No root ids need to be written to "{}"'.format(self.roots_file))
            return
        _LOG.info('Writing {} root_ids to "{}"'.format(len(root_ids), self.roots_file))
        pd = os.path.split(self.roots_file)[0]
        assure_dir_exists(pd)
        with codecs.open(self.roots_file, 'w', encoding='utf-8') as outp:
            outp.write('\n'.join([str(i) for i in root_ids]))

    def append_write_roots(self, root_ids):
        if self.garbage_bin is not None:
            self.garbage_bin.append_write_roots(root_ids)
            return
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
        if self.garbage_bin is not None:
            self.garbage_bin.write_lines(header, syn_header)
            return
        _write_taxon(header, self.id_to_line, self.id_order, self.taxon_fp)
        self.write_synonyms(syn_header)

    def append_lines(self):
        if self.garbage_bin is not None:
            self.garbage_bin.append_lines()
            return
        _append_taxon(self.id_to_line, self.id_order, self.taxon_fp)
        self.append_synonyms()

    @property
    def existing_output(self):
        if os.path.exists(self.taxon_fp):
            return self.taxon_fp
        return None

    def is_same_dest_pe(self, other):
        return self.taxon_fp == other.taxon_fp

    def add_synonym(self, el_id, line):
        if self.syn_fp is None:
            self.add_taxon(el_id, None, line)
        elif el_id:
            self.syn_by_id[el_id] = line
            self.syn_id_order.append(el_id)
        else:
            self.id_less_syn.append(line)

    def write_synonyms(self, header):
        if self.garbage_bin is not None:
            return self.garbage_bin.write_synonyms(header)
        if self.syn_fp is None:
            return
        if self.id_less_syn:
            assert not self.syn_by_id
            _write_taxon_list(header, self.id_less_syn, self.syn_fp)
        else:
            _write_taxon(header, self.syn_by_id, self.syn_id_order, self.syn_fp)

    def append_synonyms(self):
        if self.garbage_bin is not None:
            return self.garbage_bin.append_synonyms()
        if self.syn_fp is None:
            return
        if self.id_less_syn:
            assert not self.syn_by_id
            _append_taxon_list(self.id_less_syn, self.syn_fp)
        else:
            _append_taxon(self.syn_by_id, self.syn_id_order, self.syn_fp)

    def _finish_partition_from_dict(self):
        while True:
            par_id_matched = [p for p in self.id_to_child_list.keys() if p in self.id_to_el]
            if not par_id_matched:
                break
            for p in par_id_matched:
                match_el = self.id_to_el[p]
                self._transfer_children(p, match_el)
                del self.id_to_child_list[p]
        misc_el = self.garbage_bin
        for par_id in self.id_to_child_list.keys():
            self._transfer_children(par_id, misc_el)
        self.id_to_child_list = None
        for par_less_id, line in self.id_to_el:
            misc_el.add_taxon(par_less_id, None, line)
        self.id_to_line = None

    def _transfer_children(self, par_id, part_element):
        id_list = self.id_to_child_list[par_id]
        for child_id in id_list:
            self.id_to_el[child_id] = part_element
            part_element.add_taxon(child_id, par_id, self.id_to_line[child_id])
            del self.id_to_line[child_id]

    def _register_synonyms(self):
        for accept_id, i_l_list in self.syn_by_id.items():
            match_el = self.id_to_el.get(accept_id)
            if match_el is None:
                match_el = self.garbage_bin
            for col_id, line in i_l_list:
                match_el.add_synonym(col_id, line)
        self.syn_by_id = None

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

    def read_taxon_line(self, uid, par_id, line):
        if uid in self.roots_set:
            # _LOG.info("{} not in {}".format(uid, roots_set))
            match_l = [i[1] for i in self.by_roots if uid in i[0]]
            assert len(match_l) == 1
            match_el = match_l[0]
            self.id_to_el[uid] = match_el
            match_el.add_taxon(uid, par_id, line)
            if self.garbage_bin is not None:
                self.garbage_bin.add_taxon(uid, par_id, line)
        else:
            # _LOG.info("{} not in {}".format(uid, roots_set))
            if par_id:
                try:
                    par_id = int(par_id)
                except:
                    pass
            match_el = self.id_to_el.get(par_id)
            if match_el is not None:
                self.id_to_el[uid] = match_el
                match_el.add_taxon(uid, par_id, line)
            else:
                self.id_to_child_list.setdefault(par_id, []).append(uid)
                self.id_to_line[uid] = line


def create_partition_element(path_pref=None,
                             fragment=None,
                             path_suffix=None,
                             roots=None,
                             syn_filename=None,
                             taxon_filepath=None):
    from taxalotl.partitions import INP_TAXONOMY_DIRNAME
    if taxon_filepath is None:
        taxon_filepath = os.path.join(path_pref,
                                      fragment,
                                      INP_TAXONOMY_DIRNAME,
                                      path_suffix)
    if syn_filename is not None:
        syn_fp = os.path.join(os.path.split(taxon_filepath)[0], syn_filename)
    else:
        syn_fp = None
    return TaxonPartition(roots=roots,
                          taxon_fp=taxon_filepath,
                          syn_fp=syn_fp)
