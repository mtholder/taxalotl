#!/usr/bin/env python
from __future__ import print_function

import sys
from enum import IntEnum, IntFlag
from typing import Dict, List

from peyotl import (get_logger)

from .config import TaxalotlConfig
from .tree import TaxonTree
from .partitions import (PART_NAMES)
from .resource_wrapper import TaxonomyWrapper
from .taxonomic_ranks import GENUS_RANK_TO_SORTING_NUMBER, SPECIES_SORTING_NUMBER
from .tax_partition import IGNORE_SYN_TYPES

_LOG = get_logger(__name__)
out_stream = sys.stdout

SEP_NAMES = '__separator_names__.json'
SEP_MAPPING = '__separator_names_to_dir__.json'


def analyze_update_to_resources(taxalotl_config: TaxalotlConfig,
                                prev: TaxonomyWrapper,
                                curr: TaxonomyWrapper,
                                level_list: List[str]):
    m = 'Could not analyze taxonomy update because {} has not been partitioned.'
    for el in [prev, curr]:
        if not el.has_been_partitioned():
            raise RuntimeError(m.format(el.id))
    if level_list == [None]:
        level_list = PART_NAMES
    for part_name in level_list:
        analyze_update_for_level(taxalotl_config, prev, curr, part_name)


class UpdateStatus(IntFlag):
    UNCHANGED = 0
    PAR_CHANGED = 1
    NAME_CHANGED = 2
    NEW_TERMINAL = 4  # really just NEW if not combined with TERMINAL or SYNONYM
    INTERNAL = 8
    UNDIAGNOSED_CHANGE = 16
    DELETED_TERMINAL = 32  # really just DELETED if not combined with TERMINAL or SYNONYM
    SYN_DELETED = 64
    PROMOTED_FROM_SYNONYM = 128
    SUNK_TO_SYNONYM = 256
    # End flags. Start of unions
    NAME_AND_PAR_CHANGED = PAR_CHANGED | NAME_CHANGED
    NEW_INTERNAL = NEW_TERMINAL | INTERNAL
    DELETED_INTERNAL = DELETED_TERMINAL | INTERNAL
    SYNONYM = SYN_DELETED
    SYN_ADDED = SYNONYM | NEW_TERMINAL
    SYN_CHANGED = SYNONYM | UNDIAGNOSED_CHANGE
    SYN_NAME_CHANGED = SYNONYM | NAME_CHANGED
    TERMINAL_PROMOTED_FROM_SYNONYM = PROMOTED_FROM_SYNONYM | NEW_TERMINAL
    INTERNAL_PROMOTED_FROM_SYNONYM = PROMOTED_FROM_SYNONYM | NEW_INTERNAL
    TERMINAL_SUNK_TO_SYNONYM = SUNK_TO_SYNONYM | DELETED_TERMINAL
    INTERNAL_SUNK_TO_SYNONYM = SUNK_TO_SYNONYM | DELETED_INTERNAL


def del_add_set_diff(old_set, new_set):
    return old_set.difference(new_set), new_set.difference(old_set)


def del_mod_add_dict_diff(oldd, newd):
    same, d, m, a = True, {}, {}, {}
    if oldd != newd:
        for k, v in oldd.items():
            if k in newd:
                nv = newd[k]
                if nv == v:
                    continue
                m[k] = nv
                same = False
            else:
                d[k] = v
                same = False
        for k, v in newd.items():
            if k not in oldd:
                a[k] = v
                same = False
    return same, (d, m, a)


def flag_update_status(f, s, stat):
    f.update_status.setdefault(stat, []).append(s)
    if s is not None:
        s.update_status.setdefault(stat, []).append(f)


def flag_synonyms_change(f, s):
    dels, adds = del_add_set_diff(s.synonyms, f.synonyms)
    real_dels = set()
    changed = set()
    for d in dels:
        n = d.name
        matching_add = None
        for a in adds:
            if a.name == n:
                matching_add = a
        if matching_add is None:
            real_dels.add(d)
        else:
            changed.add((d, matching_add))
            adds.remove(matching_add)
    for d in real_dels:
        s.update_status.setdefault(UpdateStatus.SYN_DELETED, []).append(d)
    for d in adds:
        f.update_status.setdefault(UpdateStatus.SYN_ADDED, []).append(d)
    for p in changed:
        f.update_status.setdefault(UpdateStatus.SYN_CHANGED, []).append(p)
    return real_dels, changed, adds


def _has_syn_update(nd):
    for k, nd_list in nd.update_status.items():
        if k & UpdateStatus.SYNONYM:
            return True
    return False


def _get_nonsyn_flag_and_other(nd):
    node_status_dict = nd.update_status
    node_status_flag, other_node = None, None
    for k, nd_list in node_status_dict.items():
        if k & UpdateStatus.SYNONYM:
            continue
        if node_status_flag is None:
            assert len(nd_list) == 1
            node_status_flag, other_node = k, nd_list[0]
        else:
            node_status_flag, other_node = None, None
            break
    if node_status_flag is None:
        m = 'incorrect node_status_list len = {}\n{}\n'
        raise RuntimeError(m.format(node_status_dict, nd.__dict__))
    return node_status_flag, other_node


def _old_modified_subtree_ids(init_mod_par_set, tree):
    oldest_mod_par = set()
    seen = set()
    for pid in init_mod_par_set:
        if pid in seen:
            continue
        try:
            nd = tree.id_to_taxon[pid]
        except:
            oldest_mod_par.add(pid)
            continue
        seen.add(pid)
        nsl = _get_nonsyn_flag_and_other(nd)[0]
        added = False
        while nsl != UpdateStatus.UNCHANGED:
            seen.add(nd.id)
            try:
                nd = tree.id_to_taxon[nd.par_id]
            except:
                oldest_mod_par.add(nd.id)
                added = True
                break
            nsl = _get_nonsyn_flag_and_other(nd)[0]
        if not added:
            oldest_mod_par.add(nd.id)
    return oldest_mod_par


class NodeRankCmpResult(IntEnum):
    SAME_RANK = 0
    FIRST_HIGHER = 1
    SECOND_HIGHER = 2
    OVERLAPPING = 3


def _compare_nd_ranks(first_pair, second_pair):
    fir_tree, fir_nd = first_pair
    sec_tree, sec_nd = second_pair
    fir_rr = fir_tree.node_rank_sorting_number_range(fir_nd)
    sec_rr = sec_tree.node_rank_sorting_number_range(sec_nd)
    if fir_rr == sec_rr:
        return NodeRankCmpResult.SAME_RANK, fir_rr, sec_rr
    if fir_rr[1] > sec_rr[0]:
        return NodeRankCmpResult.FIRST_HIGHER, fir_rr, sec_rr
    if sec_rr[1] > fir_rr[0]:
        return NodeRankCmpResult.SECOND_HIGHER, fir_rr, sec_rr
    return NodeRankCmpResult.OVERLAPPING, fir_rr, sec_rr


class UpdateStatusLog(object):
    def __init__(self, prev_tree=None, curr_tree=None):
        self.in_order = []
        self.by_status_code = {}  # type: Dict[UpdateStatus, List]
        self.prev_tree, self.curr_tree = None, None
        self.set_prev_curr(prev_tree, curr_tree)
        self.written_ids = set()

    def improve_status(self, old_stat, new_stat, node_list):
        to_cull = self.by_status_code.get(old_stat, [])
        for nd in node_list:
            us = nd.update_status
            us[new_stat] = us[old_stat]
            del us[old_stat]
            to_cull.remove(nd)
        self.by_status_code.setdefault(new_stat, []).extend(node_list)

    def set_prev_curr(self, prev_tree, curr_tree):
        self.prev_tree, self.curr_tree = prev_tree, curr_tree

    def add_node(self, nd):
        node_status_flag = _get_nonsyn_flag_and_other(nd)[0]
        self.in_order.append(nd)
        self.by_status_code.setdefault(node_status_flag, []).append(nd)

    def _par_is_genus_report(self, nd, is_in_curr, tree, status, other_nd):
        assert is_in_curr
        other_tree = self.prev_tree if is_in_curr else self.curr_tree
        assert other_nd
        other_rank_range = other_tree.node_rank_sorting_number_range(other_nd)
        assert other_rank_range == (GENUS_RANK_TO_SORTING_NUMBER, GENUS_RANK_TO_SORTING_NUMBER)
        nd_c = nd.child_id_dict()
        other_c = other_nd.child_id_dict()
        common_unchanged_dict = {}
        common_mod_dict = {}
        new_child_dict = {}
        for nd_c_id, nd_child in nd_c.items():
            other_nd_child = other_c.get(nd_c_id)
            if other_nd_child is not None:
                nd_c_stat = _get_nonsyn_flag_and_other(nd_child)[0]
                if nd_c_stat == UpdateStatus.UNCHANGED:
                    common_unchanged_dict[nd_c_id] = nd_child
                else:
                    assert nd_child.id == other_nd_child.id
                    common_mod_dict[nd_c_id] = nd_child
            else:
                paired_nd = _get_nonsyn_flag_and_other(nd_child)[1]
                if paired_nd is None:
                    new_child_dict[nd_c_id] = nd_child
                else:
                    nd_c_is_spec_typed = tree.node_is_specimen_typed(nd_child)
                    if nd_c_is_spec_typed:
                        assert other_tree.node_is_specimen_typrugoed(paired_nd)
                        rank_cmp, frr, scc = _compare_nd_ranks((tree, nd_child),
                                                               (other_tree, paired_nd))
                        if rank_cmp == NodeRankCmpResult.FIRST_HIGHER:
                            if frr == (SPECIES_SORTING_NUMBER, SPECIES_SORTING_NUMBER):
                                raise NotImplementedError('hi')

        del_child_dict = {k: other_c[k] for k in other_c.keys() if k not in nd_c}
        num_unmod = len(common_unchanged_dict)
        num_mod = len(common_mod_dict)
        num_new = len(new_child_dict)
        num_del = len(del_child_dict)
        m = '{}/{}/{}/{} Retained/Modified/Added/Deleted children in GENUS {}\n'
        m = m.format(num_unmod, num_mod, num_new, num_del, str(nd))
        out_stream.write(m)

        if common_mod_dict:
            out_stream.write('  modified:\n')
            for ni, cnd in common_mod_dict.items():
                self._write_nd(cnd, indent='    ')
        if new_child_dict:
            out_stream.write('  added:\n')
            for ni, cnd in new_child_dict.items():
                self._write_nd(cnd, indent='    ')
        if del_child_dict:
            out_stream.write('  deleted:\n')
            for ni, cnd in del_child_dict.items():
                self._write_nd(cnd, indent='    ')

    def report_on_altered_contiguous_des(self, nd, is_in_curr_tree):
        assert is_in_curr_tree
        tree = self.curr_tree if is_in_curr_tree else self.prev_tree
        status, other_nd = _get_nonsyn_flag_and_other(nd)
        # None if a new node, but we aren't expecting the oldest node's
        #   parent to be new
        assert other_nd is not None
        if tree.node_is_specimen_typed(nd):
            _LOG.warn('SPEC  ')
            # self._write_nd(nd, True)
        else:
            max_sn, min_csn = tree.node_rank_sorting_number_range(nd)
            if max_sn == min_csn and max_sn == GENUS_RANK_TO_SORTING_NUMBER:
                self._par_is_genus_report(nd, is_in_curr_tree, tree, status, other_nd)
            else:
                _LOG.warn('CLADE ')
                # self._write_nd(nd, True)

    def flush(self):
        curr_tree_par_ids = set()
        prev_tree_par_ids = set()
        for status_code, node_list in self.by_status_code.items():
            if status_code == UpdateStatus.UNCHANGED:
                continue
            if status_code in [UpdateStatus.DELETED_TERMINAL, UpdateStatus.DELETED_INTERNAL]:
                target = prev_tree_par_ids
            else:
                target = curr_tree_par_ids
            for nd in node_list:
                target.add(nd.par_id)

        curr_deepest_mod_id = _old_modified_subtree_ids(curr_tree_par_ids, self.curr_tree)
        prev_deepest_mod_id = _old_modified_subtree_ids(prev_tree_par_ids, self.prev_tree)

        # emitted = set()
        # for par_id in curr_deepest_mod_id:
        #     par_nd = self.curr_tree.id_to_taxon[par_id]
        #     self.report_on_altered_contiguous_des(par_nd, True)

        status_keys = [(i.value, i) for i in self.by_status_code.keys()]
        status_keys.sort()
        status_keys = [i[1] for i in status_keys]
        status_keys.remove(UpdateStatus.TERMINAL_SUNK_TO_SYNONYM)
        status_keys.remove(UpdateStatus.INTERNAL_SUNK_TO_SYNONYM)
        for k in status_keys:
            for nd in self.by_status_code[k]:
                self._write_nd(nd)

        # Reinitialize...
        self.__init__(None, None)

    def _write_nd(self, nd, even_unchanged=False, indent=''):
        node_status_flag, other_node = _get_nonsyn_flag_and_other(nd)
        changed = node_status_flag != UpdateStatus.UNCHANGED
        nsyn_c = _has_syn_update(nd)
        osyn_c = (other_node is not None) and _has_syn_update(other_node)
        has_sunken = hasattr(nd, 'sunken')
        has_prev_container = hasattr(nd, 'prev_container')
        if not changed:
            if has_sunken or has_prev_container:
                changed = True
            if osyn_c or nsyn_c:
                changed = True
        if (not changed) and (not even_unchanged):
            return
        if nd.id in self.written_ids:
            return
        self.written_ids.add(nd.id)
        msg_template = '{}{}: {}. Previous: {}\n'
        m = msg_template.format(indent, node_status_flag.name, nd, other_node)
        out_stream.write(m)
        if has_prev_container:
            for x in nd.prev_container:
                xnsf, xon = _get_nonsyn_flag_and_other(x)
                m = msg_template.format('   prev. container: ' + indent, xnsf.name, x, xon)
                out_stream.write(m)
        if has_sunken:
            for x in nd.sunken:
                xnsf, xon = _get_nonsyn_flag_and_other(x)
                m = msg_template.format('   ' + indent, xnsf.name, x, xon)
                out_stream.write(m)
        if nsyn_c:
            self._write_syn_for_nd(nd, indent)
        if node_status_flag == UpdateStatus.UNCHANGED and osyn_c:
            self._write_syn_for_nd(other_node, indent)

    def _write_syn_for_nd(self, nd, indent):
        for f in [UpdateStatus.SYN_DELETED, UpdateStatus.SYN_ADDED, UpdateStatus.SYN_CHANGED]:
            sd = nd.update_status.get(f)
            if sd is not None:
                m = '   {}{}: {}\n'.format(indent, f.name, sd)
                out_stream.write(m)


def analyze_update_for_level(taxalotl_config: TaxalotlConfig,
                             prev: TaxonomyWrapper,
                             curr: TaxonomyWrapper,
                             part_name: str,
                             update_log: UpdateStatusLog = None):
    if update_log is None:
        update_log = UpdateStatusLog()
    fragment = taxalotl_config.get_fragment_from_part_name(part_name)
    _LOG.info('analyze_update_to_resources for {} for {} -> {}'.format(fragment, prev.id, curr.id))
    pf = prev.get_taxon_forest_for_partition(part_name)
    cf = curr.get_taxon_forest_for_partition(part_name)
    if len(pf.trees) != 1 or len(cf.trees) != 1:
        raise NotImplementedError('Analysis of multi-tree forests not supported, yet')
    prev_tree = pf.trees[0]  # type: TaxonTree
    # prev_tree.write_rank_indented(out_stream)
    curr_tree = cf.trees[0]  # type: TaxonTree
    update_log.set_prev_curr(prev_tree, curr_tree)
    prev_tree.add_num_tips_below()
    curr_tree.add_num_tips_below()
    prev_tree.add_update_fields()
    curr_tree.add_update_fields()
    print(prev_tree.root.num_tips_below, '->', curr_tree.root.num_tips_below)
    prev_syn = prev.get_parsed_synonyms_by_id(part_name, ignored_syn_types=IGNORE_SYN_TYPES)
    prev_tree.attach_parsed_synonyms_set(prev_syn)
    curr_syn = curr.get_parsed_synonyms_by_id(part_name, ignored_syn_types=IGNORE_SYN_TYPES)
    curr_tree.attach_parsed_synonyms_set(curr_syn)

    # same_synonyms, syn_dma_dicts = del_mod_add_dict_diff(prev_syn, curr_syn)
    # if same_synonyms:
    #     _LOG.info('No changes to synonyms')
    # else:
    #     d, m, a = syn_dma_dicts
    #     print(len(d), len(m), len(a))
    #     for k, v in d.items():
    #         print('Deleted SYNONYM set: {} => {}'.format(k, v))
    #     for k, v in a.items():
    #         print('Added   SYNONYM set: {} => {}'.format(k, v))
    #     for k in m.keys():
    #         ov = prev_syn[k]
    #         nv = curr_syn[k]
    #         dels, adds = del_add_set_diff(ov, nv)
    #         for d in dels:
    #             print('Deleted SYNONYM    : {} => {}'.format(k, d))
    #         for d in adds:
    #             print('Added   SYNONYM    : {} => {}'.format(k, d))
    # sys.exit('hi')
    for leaf in curr_tree.postorder():
        lid = leaf.id
        prev_nd = prev_tree.id_to_taxon.get(lid)
        if prev_nd is None:
            s = UpdateStatus.NEW_INTERNAL if leaf.children_refs else UpdateStatus.NEW_TERMINAL
            flag_update_status(leaf, None, s)
        else:
            psd = prev_nd.to_serializable_dict()
            lsd = leaf.to_serializable_dict()
            same, dma_dicts = del_mod_add_dict_diff(psd, lsd)
            if same:
                flag_update_status(leaf, prev_nd, UpdateStatus.UNCHANGED)
            else:
                d, m, a = dma_dicts
                if (not d) and (not a):
                    if len(m) == 1:
                        if 'par_id' in m:
                            flag_update_status(leaf, prev_nd, UpdateStatus.PAR_CHANGED)
                        elif 'name' in m:
                            flag_update_status(leaf, prev_nd, UpdateStatus.NAME_CHANGED)
                    elif len(m) == 2 and 'par_id' in m and 'name' in m:
                        flag_update_status(leaf, prev_nd, UpdateStatus.NAME_AND_PAR_CHANGED)
                if not leaf.update_status:
                    flag_update_status(leaf, prev_nd, UpdateStatus.UNDIAGNOSED_CHANGE)

    for prev_nd in prev_tree.preorder():
        if not prev_nd.update_status:
            is_intern = bool(prev_nd.children_refs)
            s = UpdateStatus.DELETED_INTERNAL if is_intern else UpdateStatus.DELETED_TERMINAL
            flag_update_status(prev_nd, None, s)
            update_log.add_node(prev_nd)

    del_syn, mod_syn, add_syn = set(), set(), set()
    for nd in curr_tree.preorder():
        other = _get_nonsyn_flag_and_other(nd)[-1]
        if other is not None:
            if other.synonyms != nd.synonyms:
                d, m, a = flag_synonyms_change(nd, other)
                del_syn.update(d)
                mod_syn.update(m)
                add_syn.update(a)
        update_log.add_node(nd)

    del_syn_names, mod_syn_names, add_syn_names = {}, {}, {}
    for d in del_syn:
        del_syn_names.setdefault(d.name, []).append(d)
    # for d in mod_syn:
    #     mod_syn_names.setdefault(d.name, []).append(d)
    for d in add_syn:
        add_syn_names.setdefault(d.name, []).append(d)

    for k in [UpdateStatus.NEW_INTERNAL, UpdateStatus.NEW_TERMINAL]:
        more_specific_status = set()
        for nd in update_log.by_status_code.get(k, []):
            if nd.name in del_syn_names:
                psl = del_syn_names[nd.name]
                nd.prev_syn = psl
                for prevsyn in psl:
                    prev_container = prev_tree.id_to_taxon[prevsyn.valid_tax_id]
                    if not hasattr(nd, 'prev_container'):
                        nd.prev_container = []
                    nd.prev_container.append(prev_container)
                more_specific_status.add(nd)
        ns = UpdateStatus.PROMOTED_FROM_SYNONYM | k
        update_log.improve_status(k, ns, more_specific_status)

    for k in [UpdateStatus.DELETED_TERMINAL, UpdateStatus.DELETED_INTERNAL]:
        more_specific_status = set()
        for nd in update_log.by_status_code.get(k, []):
            if nd.name in add_syn_names:
                nsl = add_syn_names[nd.name]
                for new_syn in nsl:
                    target = curr_tree.id_to_taxon[new_syn.valid_tax_id]
                    if not hasattr(target, 'sunken'):
                        target.sunken = []
                    target.sunken.append(nd)
                nd.new_syn = nsl
                more_specific_status.add(nd)
        ns = UpdateStatus.SUNK_TO_SYNONYM | k
        update_log.improve_status(k, ns, more_specific_status)

    x = update_log.flush()
    print('done')
    return x
