#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from peyotl import (add_or_append_to_dict, assure_dir_exists,
                    get_logger,
                    write_as_json)
import tempfile
import codecs
import csv
import os

_LOG = get_logger(__name__)

INP_OTT_TAXONOMY_HEADER = "uid\t|\tparent_uid\t|\tname\t|\trank\t|\t\n"
INP_FLAGGED_OTT_TAXONOMY_HEADER = "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tflags\t|\t\n"
INP_OTT_SYNONYMS_HEADER = "uid\t|\tname\t|\ttype\t|\t\n"


def write_ott_taxonomy_tsv(out_fp,
                           root_nodes,
                           id_to_par,
                           id_to_children,
                           id_to_rank,
                           id_to_name,
                           has_syn_dict,
                           details_log,
                           extinct_known=None):
    """If has_syn_dict is provided, then a list of the IDs that occur in that dict
    is returned in the order that the IDs were written to taxonomy file. This 
    allows for the synonyms.tsv file to be written in a similar order, which makes browsing it
    easier.
    """
    if has_syn_dict is None:
        has_syn_dict = {}
    syn_id_order = []
    num_tips_written = 0
    num_internals_written = 0
    rn = list(root_nodes)
    rn.sort()
    header = INP_OTT_TAXONOMY_HEADER if extinct_known is None else INP_FLAGGED_OTT_TAXONOMY_HEADER
    with codecs.open(out_fp, 'w', encoding='utf-8') as out:
        out.write(header)
        # need to print id, parent id, and name
        for root_id in rn:
            stack = [root_id]
            while stack:
                curr_id = stack.pop()
                try:
                    if curr_id in has_syn_dict:
                        syn_id_order.append(curr_id)
                    name = id_to_name[curr_id]
                    par_id = id_to_par.get(curr_id)
                    if par_id is None:
                        spar_id = ''
                    else:
                        spar_id = str(par_id)
                    rank = id_to_rank.get(curr_id, '')
                    children = id_to_children.get(curr_id)
                    if children:
                        num_internals_written += 1
                        stack.extend(children)
                    else:
                        num_tips_written += 1
                    fields = [str(curr_id), spar_id, name, rank, '']
                    if extinct_known is not None:
                        ev = extinct_known.get(curr_id)
                        if ev:
                            fields[-1] = 'extinct'
                        fields.append('')
                    try:
                        out.write(u'{}\n'.format(u'\t|\t'.join(fields)))
                    except:
                        _LOG.exception("error serializing {}".format(repr(fields)))
                except:
                    _LOG.error("Error writing taxon_id {}".format(curr_id))
                    raise
    details_log['num_tips_written'] = num_tips_written
    details_log['num_internals_written'] = num_internals_written
    return syn_id_order


def write_ott_synonyms_tsv(out_fp,
                           id_to_name_name_type_list,
                           id_order,
                           details_log):
    num_syn_written = 0
    with codecs.open(out_fp, 'w', encoding='utf-8') as out:
        out.write(INP_OTT_SYNONYMS_HEADER)
        for nd_id in id_order:
            syn_list = id_to_name_name_type_list[nd_id]
            for name, name_type, syn_id in syn_list:
                num_syn_written += 1
                out.write(u'{}\n'.format('\t|\t'.join([str(nd_id), name, name_type, ''])))
    details_log['num_synonyms_written'] = num_syn_written
    details_log['num_ids_with_synonyms_written'] = len(id_order)


def write_ott_forwards(out_fp, forwarded_dict):
    with codecs.open(out_fp, 'w', encoding='utf-8') as out:
        for key, value in forwarded_dict.items():
            out.write('{}\t{}\n'.format(key, value))


def write_ncbi_details_json(fp, details_log):
    write_as_json(details_log, fp, indent=2)


def read_taxonomy_to_get_id_to_name(tax_dir, id_coercion=int):
    ncbi_to_name = {}
    i = 0
    fp = os.path.join(tax_dir, 'taxonomy.tsv')
    with codecs.open(fp, 'r', encoding='utf-8') as inp:
        reader = csv.reader(inp, delimiter='\t')
        header = reader.next()
        uidx = header.index('uid')
        namex = header.index('name')
        for row in reader:
            uid = id_coercion(row[uidx])
            name = row[namex]
            if name is not None:
                ncbi_to_name[uid] = name
                i += 1
                if i % 200000 == 0:
                    _LOG.info("{} {} {}".format(i, uid, name))
    return ncbi_to_name


class OTTTaxon(object):
    def __init__(self, interim_taxonomy_format_line, line_num):
        self.line_num = line_num
        line = interim_taxonomy_format_line
        try:
            ls = line.split('\t|\t')
            assert len(ls) == 8 and ls[-1] == '\n'
        except:
            _LOG.exception("Error reading line {}:\n{}".format(line_num, line))
            raise
        self.id = int(ls[0])
        if ls[1]:
            self.par_id = int(ls[1])
        else:
            self.par_id = None
        self.name, self.rank, self.uniqname = ls[2], ls[3], ls[5]
        self.flags = set(ls[6].split(','))
        sel = ls[4].split(',')
        d = {}
        for el in sel:
            src, sid = el.split(':')
            try:
                sid = int(sid)
            except:
                pass
            d.setdefault(src, set()).add(sid)
        self.src_dict = d

    @property
    def name_that_is_unique(self):
        return self.uniqname if self.uniqname else self.name


def read_taxonomy_to_get_id_to_fields(tax_dir):
    fp = os.path.join(tax_dir, 'taxonomy.tsv')
    fields = ['uid', 'parent_uid', 'name', 'rank', 'sourceinfo', 'uniqname', 'flags', '\n']
    expected_header = '\t|\t'.join(fields)
    if not os.path.exists(fp):
        return {}
    with codecs.open(fp, 'r', encoding='utf-8') as inp:
        iinp = iter(inp)
        header = iinp.next()
        assert header == expected_header
        id_to_obj = {}
        for n, line in enumerate(iinp):
            obj = OTTTaxon(line, line_num=1 + n)
            oid = obj.id
            assert oid not in id_to_obj
            id_to_obj[oid] = obj
        return id_to_obj


def read_taxonomy_to_get_single_taxon(tax_dir, root_id):
    sri = str(root_id)
    fp = os.path.join(tax_dir, 'taxonomy.tsv')
    fields = ['uid', 'parent_uid', 'name', 'rank', 'sourceinfo', 'uniqname', 'flags', '\n']
    expected_header = '\t|\t'.join(fields)
    with codecs.open(fp, 'r', encoding='utf-8') as inp:
        iinp = iter(inp)
        header = iinp.next()
        assert header == expected_header
        for n, line in enumerate(iinp):
            if not line.startswith(sri):
                continue
            obj = OTTTaxon(line, line_num=1 + n)
            if root_id == obj.id:
                return obj


class InterimTaxonomyData(object):
    def __init__(self):
        self.about = {}
        self.details_log = {}
        self.forwards = {}  # from old ID to new ID
        self.to_par = {}  # ID -> parent ID or None
        self.to_children = {}  # ID to list of children IDs
        self.to_rank = {}  # ID -> rank string
        self.root_nodes = set()  # set of IDs
        self.to_name = {}  # ID -> name
        self.name_to_ids = {}  # name to
        self.synonyms = {}
        self.repeated_names = set()
        self.extinct_known = None
        self.syn_id_to_valid = None
        self.extra_blob = None

    def finalize(self):
        self.details_log['num_forwards'] = len(self.forwards)
        self.details_log['num_nodes'] = len(self.to_par)
        self.details_log['num_distinct_names'] = len(self.name_to_ids)
        self.details_log['num_ids_with_synonyms'] = len(self.synonyms)

    def register_id_and_name(self, taxon_id, name):
        self.to_name[taxon_id] = name
        if add_or_append_to_dict(self.name_to_ids, name, taxon_id):
            self.repeated_names.add(name)

    def register_synonym(self, valid_id, syn_name, name_type, syn_id=None):
        assert valid_id != syn_id
        self.synonyms.setdefault(valid_id, []).append((syn_name, name_type, syn_id))

    def fix_synonym(self, valid_id, old_valid, syn_id):
        assert valid_id != old_valid
        ol = self.synonyms.get(old_valid, [])
        matching = []
        for el in ol:
            if el[2] == syn_id:
                matching.append(el)
        nl = self.synonyms.setdefault(valid_id, [])
        for el in matching:
            ol.remove(el)
            nl.append(el)

    def write_ott_taxonomy_tsv(self, fp):
        return write_ott_taxonomy_tsv(fp,
                                      self.root_nodes,
                                      self.to_par,
                                      self.to_children,
                                      self.to_rank,
                                      self.to_name,
                                      self.synonyms,
                                      self.details_log,
                                      self.extinct_known)

    def write_to_dir(self, destination):
        # Write out in OTT form
        d = tempfile.mkdtemp()
        fn = ['taxonomy.tsv',
              'synonyms.tsv',
              'forwards.tsv',
              'about.json',
              'details.json']
        try:
            syn_order = self.write_ott_taxonomy_tsv(os.path.join(d, 'taxonomy.tsv'))
            write_ott_synonyms_tsv(os.path.join(d, 'synonyms.tsv'),
                                   self.synonyms,
                                   syn_order,
                                   self.details_log)
            if self.forwards:
                write_ott_forwards(os.path.join(d, 'forwards.tsv'), self.forwards)

            about_fp = os.path.join(d, 'about.json')
            write_as_json(self.about, about_fp, indent=2)
            self.finalize()
            write_ncbi_details_json(os.path.join(d, 'details.json'),
                                    self.details_log)
        except:
            for f in fn:
                tf = os.path.join(d, f)
                if os.path.exists(tf):
                    try:
                        os.remove(tf)
                    except:
                        pass
            try:
                os.rmdir(d)
            except:
                pass
            raise
        assure_dir_exists(destination)
        for f in fn:
            sfp = os.path.join(d, f)
            if os.path.exists(sfp):
                dfp = os.path.join(destination, f)
                os.rename(sfp, dfp)
        os.rmdir(d)

    def del_ids(self, id_list):
        to_name = self.to_name
        to_par = self.to_par
        to_children = self.to_children
        for taxon_id in id_list:
            if taxon_id in to_name:
                del to_name[taxon_id]
            if taxon_id in to_children:
                del to_children[taxon_id]
            pid = to_par.get(taxon_id)
            if pid:
                del to_par[taxon_id]
                pc = to_children.get(pid)
                try:
                    if pc:
                        pc.remove(taxon_id)
                except:
                    pass
