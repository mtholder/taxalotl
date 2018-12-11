#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import csv
import io
import os
import tempfile

from peyotl import (add_or_append_to_dict, assure_dir_exists,
                    get_logger,
                    shorter_fp_form,
                    write_as_json)

_LOG = get_logger(__name__)

INP_OTT_TAXONOMY_HEADER = "uid\t|\tparent_uid\t|\tname\t|\trank\t|\t\n"
INP_FLAGGED_OTT_TAXONOMY_HEADER = "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tflags\t|\t\n"
INP_FLAGGED_OTT_TAXONOMY_NO_TRAIL_HEADER = "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tflags\n"
INP_OTT_SYNONYMS_HEADER = "uid\t|\tname\t|\ttype\t|\t\n"
FULL_OTT_HEADER = "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t\n"


def _parse_synonyms(tax_part):  # type (TaxonPartition) -> None
    syn_fp = tax_part.input_synonyms_filepath
    tax_part.syn_header = ''
    if not os.path.exists(syn_fp):
        return
    _LOG.debug('parsing synonyms from "{}" ...'.format(syn_fp))
    with io.open(syn_fp, 'rU', encoding='utf-8') as inp:
        iinp = iter(inp)
        try:
            tax_part.syn_header = next(iinp)
        except StopIteration:
            return
        shs = tax_part.syn_header.split('\t|\t')
        if shs[0] == 'uid':
            uid_ind = 0
        elif shs[1] == 'uid':
            uid_ind = 1
        else:
            raise ValueError("Expected one of the first 2 columns of an OTT formatted "
                             "synonyms file to be 'uid'. Problem reading: {}".format(syn_fp))
        for n, line in enumerate(iinp):
            ls = line.split('\t|\t')
            if n > 0 and n % 1000 == 0:
                _LOG.debug(' read synonym {:7} from "{}"'.format(n, syn_fp))
            try:
                accept_id = ls[uid_ind]
                try:
                    accept_id = int(accept_id)
                except:
                    pass
                tax_part.add_synonym(accept_id, syn_id=None, line=line)
            except:
                _LOG.exception("Exception parsing line {}:\n{}".format(1 + n, line))
                raise


def _parse_taxa(tax_part):  # type (TaxonPartition) -> None
    complete_taxon_fp = tax_part.tax_fp
    tax_part.taxon_header = ''
    if not os.path.exists(complete_taxon_fp):
        return
    ptp = shorter_fp_form(complete_taxon_fp)
    _LOG.debug('parsing taxa from "{}" ...'.format(ptp))
    with io.open(complete_taxon_fp, 'rU', encoding='utf-8') as inp:
        iinp = iter(inp)
        try:
            tax_part.taxon_header = next(iinp)
        except StopIteration:
            return
        for n, line in enumerate(iinp):
            ls = line.split('\t|\t')
            if n > 0 and n % 10000 == 0:
                _LOG.debug(' read taxon {:<7} from "{}" ...'.format(n, ptp))
            try:
                uid, par_id = ls[0], ls[1]
                try:
                    uid = int(uid)
                except:
                    pass
                tax_part.read_taxon_line(uid, par_id, line)
            except:
                _LOG.exception("Exception parsing line {}:\n{}".format(1 + n, line))
                raise


def partition_ott_by_root_id(tax_part):  # type (TaxonPartition) -> None
    _parse_synonyms(tax_part)
    _parse_taxa(tax_part)


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
    with io.open(out_fp, 'w', encoding='utf-8') as out:
        out.write(header)
        # need to print id, parent id, and name
        for root_id in rn:
            stack = [root_id]
            while stack:
                curr_id = stack.pop()
                if curr_id in has_syn_dict:
                    syn_id_order.append(curr_id)
                try:
                    name = id_to_name[curr_id]
                except KeyError:
                    _LOG.warn('Could not find a name for ID "{}"'.format(curr_id))
                    continue
                try:
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
    with io.open(out_fp, 'w', encoding='utf-8') as out:
        out.write(INP_OTT_SYNONYMS_HEADER)
        for nd_id in id_order:
            syn_list = id_to_name_name_type_list[nd_id]
            for name, name_type, syn_id in syn_list:
                num_syn_written += 1
                out.write(u'{}\n'.format('\t|\t'.join([str(nd_id), name, name_type, ''])))
    details_log['num_synonyms_written'] = num_syn_written
    details_log['num_ids_with_synonyms_written'] = len(id_order)


def write_ott_forwards(out_fp, forwarded_dict):
    with io.open(out_fp, 'w', encoding='utf-8') as out:
        for key, value in forwarded_dict.items():
            out.write('{}\t{}\n'.format(key, value))


def write_ncbi_details_json(fp, details_log):
    write_as_json(details_log, fp, indent=2)


def read_taxonomy_to_get_id_to_name(tax_dir, id_coercion=int):
    ncbi_to_name = {}
    i = 0
    fp = os.path.join(tax_dir, 'taxonomy.tsv')
    with io.open(fp, 'rU', encoding='utf-8') as inp:
        reader = csv.reader(inp, delimiter='\t')
        header = next(reader)
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


def full_ott_line_parser(taxon, line):
    try:
        ls = line.split('\t|\t')
        assert ls[-1] == '\n'
    except:
        _LOG.exception("Error reading line {}:\n{}".format(taxon.line_num, line))
        raise
    taxon.id = int(ls[0])
    if ls[1]:
        taxon.par_id = int(ls[1])
    else:
        taxon.par_id = None
    taxon.name = ls[2]
    if len(ls) <= 4:
        return
    if ls[3]:
        taxon.rank = ls[3]
    if len(ls) == 5:
        return
    if ls[4]:
        sel = ls[4].split(',')
        d = {}
        for el in sel:
            src, sid = el.split(':')
            try:
                sid = int(sid)
            except:
                pass
            d.setdefault(src, set()).add(sid)
        if d:
            taxon.src_dict = d
    if len(ls) == 6:
        return
    if ls[5]:
        taxon.uniqname = ls[5]
    if len(ls) > 7:
        if ls[6]:
            taxon.flags = set(ls[6].split(','))


def flag_after_rank_parser(taxon, line):
    try:
        ls = line.split('\t|\t')
        if len(ls) == 5:
            assert ls[4].endswith('\n')
            ls[4] = ls[4].strip()
        else:
            assert len(ls) == 6 and ls[-1] == '\n'
    except:
        _LOG.exception("Error reading line {}:\n{}".format(taxon.line_num, line))
        raise
    taxon.id = int(ls[0])
    if ls[1]:
        taxon.par_id = int(ls[1])
    else:
        taxon.par_id = None
    taxon.name = ls[2]
    if ls[3]:
        taxon.rank = ls[3]
    if ls[4]:
        taxon.flags = set(ls[4].split(','))


HEADER_TO_LINE_PARSER = {FULL_OTT_HEADER: full_ott_line_parser,
                         INP_OTT_TAXONOMY_HEADER: full_ott_line_parser,
                         INP_FLAGGED_OTT_TAXONOMY_HEADER: flag_after_rank_parser,
                         INP_FLAGGED_OTT_TAXONOMY_NO_TRAIL_HEADER: flag_after_rank_parser,
                         }


class OTTTaxon(object):
    _DATT = ('id', 'par_id', 'name', 'rank', 'src_dict', 'flags', 'uniqname')

    def __init__(self, line=None, line_num='<unknown>', line_parser=full_ott_line_parser, d=None):
        self.id, self.par_id, self.name, self.rank = None, None, None, None
        self.src_dict, self.flags, self.uniqname = None, None, None
        if d is not None:
            self.from_serializable_dict(d)
        else:
            self.line_num = line_num
            self.line = line
            line_parser(self, line)

    @property
    def name_that_is_unique(self):
        return self.uniqname if self.uniqname else self.name

    def from_serializable_dict(self, d):
        for k, v in d.items():
            if k == 'flags':
                v = set(v)
            elif k == 'src_dict':
                v = {sk: set(sv) for sk, sv in v.items()}
            setattr(self, k, v)

    def to_serializable_dict(self):
        d = {}
        for k in OTTTaxon._DATT:
            v = getattr(self, k, None)
            if v is not None:
                if k == 'id' or k == 'par_id':
                    d[k] = v
                elif v:
                    if k == 'flags':
                        if len(v):
                            d[k] = list(v)
                            d[k].sort()
                    elif k == 'src_dict':
                        ds = {}
                        for sk, ss in v.items():
                            sl = list(ss)
                            sl.sort()
                            ds[sk] = sl
                        d[k] = ds
                    else:
                        d[k] = v
        return d


def write_indented_subtree(out, node, indent_level):
    out.write('{}{} (id={})\n'.format('  ' * indent_level,
                                      node.name_that_is_unique,
                                      node.id))
    if node.children_refs:
        for c in node.children_refs:
            write_indented_subtree(out, c, indent_level=1 + indent_level)


class TaxonTree(object):
    def __init__(self, root_id, id_to_children_ids, id_to_taxon):
        self.root = id_to_taxon[root_id]
        self.root.parent_ref = None
        self.id_to_taxon = {}
        to_process = {root_id}
        while to_process:
            curr_nd_id = to_process.pop()
            curr_taxon = id_to_taxon[curr_nd_id]
            self.id_to_taxon[curr_nd_id] = curr_taxon
            curr_children_ids = id_to_children_ids.get(curr_nd_id)
            if curr_children_ids:
                curr_taxon.children_refs = [id_to_taxon[i] for i in curr_children_ids]
                for ct in curr_taxon.children_refs:
                    ct.parent_ref = curr_taxon
                to_process.update(curr_children_ids)
            else:
                curr_taxon.children_refs = None

    def get_taxon(self, uid):
        return self.id_to_taxon.get(uid)


class TaxonForest(object):
    def __init__(self, id_to_taxon):
        id_to_par = {}
        id_to_children = {}
        for taxon_id, taxon in id_to_taxon.items():
            id_to_par[taxon_id] = taxon.par_id
            id_to_children.setdefault(taxon.par_id, set()).add(taxon_id)
        root_pars = set(id_to_children.keys()) - set(id_to_par.keys())
        roots = set()
        for rp in root_pars:
            roots.update(id_to_children[rp])
        self.roots = {}
        for r in roots:
            self.roots[r] = TaxonTree(root_id=r,
                                      id_to_children_ids=id_to_children,
                                      id_to_taxon=id_to_taxon)

    def write_indented(self, out):
        for r in self.roots.values():
            write_indented_subtree(out, r.root, indent_level=0)

    @property
    def trees(self):
        return tuple(self.roots.values())

    def get_taxon(self, uid):
        for v in self.roots.values():
            t = v.get_taxon(uid)
            if t is not None:
                return t
        return None


# noinspection PyTypeChecker
def read_taxonomy_to_get_id_to_fields(tax_dir):
    fp = os.path.join(tax_dir, 'taxonomy.tsv')
    fields = ['uid', 'parent_uid', 'name', 'rank', 'sourceinfo', 'uniqname', 'flags', '\n']
    expected_header = '\t|\t'.join(fields)
    if not os.path.exists(fp):
        return {}
    with io.open(fp, 'rU', encoding='utf-8') as inp:
        iinp = iter(inp)
        header = next(iinp)
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
    with io.open(fp, 'rU', encoding='utf-8') as inp:
        iinp = iter(inp)
        header = next(iinp)
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
