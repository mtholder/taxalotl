#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from typing import Dict, List
from .taxon import Taxon

from .taxonomic_ranks import MINIMUM_SORTING_NUMBER, SPECIES_SORTING_NUMBER


def write_indented_subtree(out, node, indent_level):
    out.write('{}{} (id={})\n'.format('  ' * indent_level,
                                      node.name_that_is_unique,
                                      node.id))
    if node.children_refs:
        for c in node.children_refs:
            write_indented_subtree(out, c, indent_level=1 + indent_level)


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
        rank_sn_to_indent_n = {i: n for n, i in enumerate(ranks_sn_list)}

        for nd in self.preorder():
            indent_num = 2 * (
                    len(ranks_sn_list) - 1 - rank_sn_to_indent_n[nd.best_rank_sort_number])
            out_stream.write(
                '{}"{}" (id={}, rank={})\n'.format(' ' * indent_num, nd.name, nd.id, nd.rank))

    def add_best_guess_rank_sort_number(self):
        for nd in self.postorder():
            if hasattr(nd, 'best_rank_sort_number'):
                continue
            if nd.id == 1138017:
                print(nd.id)
            nrr = self.node_rank_sorting_number_range(nd)
            nd.best_rank_sort_number = nrr[1]
            if not nd.children_refs:
                continue
            for c in nd.children_refs:
                if c.best_rank_sort_number >= nd.best_rank_sort_number:
                    raise ValueError('rank conflict {} and child {}'.format(nd, c))

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
            if hr is None or csn > hr:
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
            if nd.children_refs:
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
            taxon.update_status = []

    def to_root_gen(self, taxon):
        if taxon is None:
            return
        yield taxon
        while taxon.par_id:
            if taxon is self.root:
                return
            taxon = self.id_to_taxon[taxon.par_id]
            yield taxon


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
