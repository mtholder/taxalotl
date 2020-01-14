#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from typing import Dict, List
from peyotl import get_logger
from .taxon import Taxon
from .taxonomic_ranks import (GENUS_SORTING_NUMBER,
                              MINIMUM_SORTING_NUMBER,
                              SPECIES_SORTING_NUMBER)

_LOG = get_logger(__name__)


def write_indented_subtree(out, node, indent_level):
    fmsd = node.formatted_src_dict()
    f = 'flags={}'.format(', '.join(node.sorted_flags)) if node.sorted_flags else ''
    s = 'src={}'.format(fmsd) if fmsd else ''
    m = '{}{}\t|\tid={}\t|\trank={}\t|\t{}\t|\t{}\n'
    out.write(m.format('    ' * indent_level, node.name_that_is_unique, node.id, node.rank, s, f))
    if node.children_refs:
        sortable = []
        for c in node.children_refs:
            sortable.append((c.name_that_is_unique, c))
        sortable.sort()
        for el in sortable:
            write_indented_subtree(out, el[1], indent_level=1 + indent_level)


class TaxonTree(object):
    def __init__(self,
                 root_id: int,
                 id_to_children_ids: Dict[int, List[int]],
                 id_to_taxon: Dict[int, Taxon],
                 taxon_partition=None):
        self.taxon_partition = taxon_partition
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

    def write_rank_indented(self, out_stream):
        self.add_best_guess_rank_sort_number()
        rr = self.root.best_rank_sort_number
        ranks_sn_set = {i.best_rank_sort_number for i in self.preorder()}
        ranks_sn_list = list(ranks_sn_set)
        ranks_sn_list.sort()
        rank_sn2indent_n = {i: n for n, i in enumerate(ranks_sn_list)}
        lrsl = len(ranks_sn_list)
        for n, nd in enumerate(self.preorder()):
            # if n > 20:
            #    return
            if n == 0:
                nd.child_indent = ''
                indent = ''
            else:
                par = self.get_taxon(nd.par_id)
                par_pref = par.child_indent
                is_first_child = nd is par.children_refs[0]
                is_last_child = nd is par.children_refs[-1]
                pni = len(par_pref)
                indent_num = 2 * (lrsl - 1 - rank_sn2indent_n[nd.best_rank_sort_number])
                indent_num -= (pni + 1)
                tail = '-' * indent_num
                if is_last_child:
                    my_prompt = '\\' + tail
                    cs = ' '
                else:
                    my_prompt = '+' + tail
                    cs = '|'
                nd.child_indent = '{}{}{}'.format(par_pref, cs, ' ' * indent_num)
                indent = '{}{}'.format(par_pref, my_prompt)
            out_stream.write('{}{}'.format(indent, nd.terse_descrip()))
        for nd in self.preorder():
            del nd.child_indent

    def does_first_contain_second(self, other_genus, other):
        while True:
            if other is other_genus:
                return True
            try:
                other = self.id_to_taxon[other.par_id]
            except:
                return False

    def find_genus_for_alpha(self, nd):
        if nd.best_rank_sort_number == GENUS_SORTING_NUMBER:
            return nd
        try:
            p = self.id_to_taxon[nd.par_id]
        except:
            return None
        # noinspection PyTypeChecker
        return self.find_genus_for_alpha(p)

    def add_best_guess_rank_sort_number(self):
        for nd in self.postorder():
            if hasattr(nd, 'best_rank_sort_number'):
                continue
            nrr = self.node_rank_sorting_number_range(nd)
            nd.best_rank_sort_number = nrr[1]
            if not nd.children_refs:
                continue
            for c in nd.children_refs:
                while c.best_rank_sort_number >= nd.best_rank_sort_number:
                    _LOG.warn('rank conflict {} and child {}, bumping parent up...'.format(nd, c))
                    nd.best_rank_sort_number += 1

    def _get_highest_child_rank(self, nd):
        if not nd.children_refs:
            return None
        hr = None
        for c in nd.children_refs:
            csn = c.rank_sorting_number()
            if csn is None:
                csn = self._get_highest_child_rank(c)
            if csn is None:
                continue
            if hr is None or csn + 1 > hr:
                hr = csn + 1
        return hr

    def _get_lowest_anc_rank_sorting_number(self, nd):
        try:
            par = self.id_to_taxon[nd.par_id]
        except:
            return None
        psn = par.rank_sorting_number()
        if psn is None:
            return self._get_lowest_anc_rank_sorting_number(par) - 1
        return psn

    def node_rank_sorting_number_range(self, nd):
        rsn = nd.rank_sorting_number()
        if rsn is not None:
            return rsn, rsn
        csn = self._get_highest_child_rank(nd)
        if csn is None:
            csn = MINIMUM_SORTING_NUMBER - 1
        asn = self._get_lowest_anc_rank_sorting_number(nd)
        if asn is None:
            asn = self._get_lowest_anc_rank_sorting_number(nd)
            raise ValueError("no rank or anc rank for {}".format(repr(nd)))
        return asn - 1, csn + 1

    def node_is_higher_taxon(self, nd):
        if self.node_is_specimen_typed(nd):
            return False

    def node_is_specimen_typed(self, nd):
        """Return True for taxa at or below species level."""
        rsn = nd.rank_sorting_number()
        if rsn is not None:
            return rsn <= SPECIES_SORTING_NUMBER
        csn = self._get_highest_child_rank(nd)
        if csn is None:
            raise ValueError("no rank or des rank for {}".format(repr(nd)))
        if csn >= SPECIES_SORTING_NUMBER:
            return False
        asn = self._get_lowest_anc_rank_sorting_number(nd)
        if asn is None:
            raise ValueError("no rank or anc rank for {}".format(repr(nd)))
        if asn <= SPECIES_SORTING_NUMBER:
            return True
        raise ValueError("Could not narrow rank for {}".format(repr(nd)))

    def get_taxon(self, uid):
        return self.id_to_taxon.get(uid)

    def postorder(self) -> Taxon:
        for t in reversed(list(self.preorder())):
            yield t

    def leaves(self) -> Taxon:
        for nd in self.preorder():
            if not nd.children_refs:
                yield nd

    def preorder(self) -> Taxon:
        curr = self.root
        if not curr:
            return
        yield curr
        to_process = []
        if curr.children_refs:
            to_process.append(list(curr.children_refs))
        while to_process:
            cl = to_process[-1]
            if not cl:
                to_process.pop(-1)
            else:
                curr = cl.pop(0)
                if curr.children_refs:
                    to_process.append(list(curr.children_refs))
                yield curr

    def add_num_tips_below(self):
        for taxon in self.postorder():
            if taxon.children_refs is None:
                taxon.num_tips_below = 1
            else:
                taxon.num_tips_below = sum([i.num_tips_below for i in taxon.children_refs])

    def add_update_fields(self):
        for taxon in self.preorder():
            taxon.update_status = {}

    def to_root_gen(self, taxon):
        if taxon is None:
            return
        yield taxon
        while taxon.par_id:
            if taxon is self.root:
                return
            taxon = self.id_to_taxon[taxon.par_id]
            yield taxon

    def attach_parsed_synonyms_set(self, syn_dict, warn_missing_target=True):
        for uid, synonym_set in syn_dict.items():
            try:
                taxon = self.id_to_taxon[uid]
            except:
                if warn_missing_target:
                    _LOG.warn('could not find target taxon for {}'.format(synonym_set))
                continue
            taxon.synonyms = synonym_set


class TaxonForest(object):
    def __init__(self, id_to_taxon, taxon_partition=None):
        self.taxon_partition = taxon_partition
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
                                      id_to_taxon=id_to_taxon,
                                      taxon_partition=taxon_partition)

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
