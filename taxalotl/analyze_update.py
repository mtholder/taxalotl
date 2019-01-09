#!/usr/bin/env python
from __future__ import print_function

import sys
from enum import Enum
from typing import Dict, List

from peyotl import (get_logger)

from taxalotl.config import TaxalotlConfig
from taxalotl.ott_schema import TaxonTree
from taxalotl.partitions import (PART_NAMES)
from taxalotl.resource_wrapper import TaxonomyWrapper
from taxalotl.ott_schema import _RANK_TO_SORTING_NUMBER

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


class UpdateStatus(Enum):
    UNCHANGED = 0
    PAR_CHANGED = 1
    NAME_CHANGED = 2
    NAME_AND_PAR_CHANGED = 3
    NEW_TERMINAL = 4
    NEW_ND_FLAG = 4
    INTERNAL_ND = 8
    NEW_INTERNAL = 12
    UNDIAGNOSED_CHANGE = 16
    DELETED_ND_FLAG = 32
    DELETED_TERMINAL = 32
    DELETED_INTERNAL = 40


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
    f.update_status.append([stat, s])
    if s is not None:
        s.update_status.append([stat, f])

def _get_flag_and_other(nd):
    node_status_list = nd.update_status
    if len(node_status_list) != 1:
        m = 'incorrect node_status_list len = {}\n{}\n'
        raise RuntimeError(m.format(node_status_list, nd.__dict__))
    node_status = node_status_list[0]
    node_status_flag = node_status[0]
    other_node = node_status[1]
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
        nsl = _get_flag_and_other(nd)[0]
        added = False
        while nsl != UpdateStatus.UNCHANGED:
            seen.add(nd.id)
            try:
                nd = tree.id_to_taxon[nd.par_id]
            except:
                oldest_mod_par.add(nd.id)
                added = True
                break
            nsl = _get_flag_and_other(nd)[0]
        if not added:
            oldest_mod_par.add(nd.id)
    return oldest_mod_par

class UpdateStatusLog(object):
    def __init__(self, prev_tree=None, curr_tree=None):
        self.in_order = []
        self.by_status_code = {} # type: Dict[UpdateStatus, List]
        self.prev_tree, self.curr_tree = None, None
        self.set_prev_curr(prev_tree, curr_tree)

    def set_prev_curr(self, prev_tree, curr_tree):
        self.prev_tree, self.curr_tree = prev_tree, curr_tree

    def add_node(self, nd):
        node_status_flag = _get_flag_and_other(nd)[0]
        self.in_order.append(nd)
        self.by_status_code.setdefault(node_status_flag, []).append(nd)

    def report_on_altered_contiguous_des(self, nd, is_in_curr_tree):
        assert is_in_curr_tree
        tree = self.curr_tree if is_in_curr_tree else self.prev_tree
        status, other_nd = _get_flag_and_other(nd)
        # None if a new node, but we aren't expecting the oldest node's
        #   parent to be new
        assert other_nd is not None
        curr_written = set()
        prev_written = set()
        if tree.node_is_specimen_typed(nd):
            out_stream.write('SPEC  ')
            self._write_nd(nd, True)
            curr_written.add(nd.id)
        else:
            max_sn, min_csn = tree.node_rank_sorting_number_range(nd)
            if max_sn == min_csn and max_sn == _RANK_TO_SORTING_NUMBER['genus']:
                out_stream.write('GENUS ')
            else:
                out_stream.write('CLADE ')

            self._write_nd(nd, True)
            curr_written.add(nd.id)

        return curr_written, prev_written

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

        emitted = set()
        for par_id in curr_deepest_mod_id:
            par_nd = self.curr_tree.id_to_taxon[par_id]
            curr_e, prev_e = self.report_on_altered_contiguous_des(par_nd, True)
            emitted.update(curr_e)
            emitted.update(prev_e)


        #     self._write_nd(nd)
        status_keys = [(i.value, i) for i in self.by_status_code.keys()]
        status_keys.sort()
        status_keys = [i[1] for i in status_keys]
        for k in status_keys:
            for nd in self.by_status_code[k]:
                self._write_nd(nd)
        # Reinitialize...
        self.__init__(None, None)

    def _write_nd(self, nd, even_unchanged=False):
        node_status_flag, other_node = _get_flag_and_other(nd)
        if node_status_flag == UpdateStatus.UNCHANGED and not even_unchanged:
            return
        m = '{}: {}. Previous: {}\n'.format(node_status_flag, nd, other_node)
        out_stream.write(m)


def analyze_update_for_level(taxalotl_config: TaxalotlConfig,
                             prev: TaxonomyWrapper,
                             curr: TaxonomyWrapper,
                             part_name: str,
                             update_log: UpdateStatusLog=None):
    if update_log is None:
        update_log = UpdateStatusLog()
    fragment = taxalotl_config.get_fragment_from_part_name(part_name)
    print('analyze_update_to_resources for {}'.format(fragment))
    pf = prev.get_taxon_forest_for_partition(part_name)
    cf = curr.get_taxon_forest_for_partition(part_name)
    if len(pf.trees) != 1 or len(cf.trees) != 1:
        raise NotImplementedError('Analysis of multi-tree forests not supported, yet')
    prev_tree = pf.trees[0]  # type: TaxonTree
    curr_tree = cf.trees[0]  # type: TaxonTree
    update_log.set_prev_curr(prev_tree, curr_tree)
    prev_tree.add_num_tips_below()
    curr_tree.add_num_tips_below()
    prev_tree.add_update_fields()
    curr_tree.add_update_fields()
    print(prev_tree.root.num_tips_below, '->', curr_tree.root.num_tips_below)
    prev_syn = prev_tree.taxon_partition.synonyms_by_id
    curr_syn = curr_tree.taxon_partition.synonyms_by_id
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
    for nd in curr_tree.preorder():
        update_log.add_node(nd)

    for prev_nd in prev_tree.preorder():
        if not prev_nd.update_status:
            is_intern = bool(prev_nd.children_refs)
            s = UpdateStatus.DELETED_INTERNAL if is_intern else UpdateStatus.DELETED_TERMINAL
            flag_update_status(prev_nd, None, s)
            update_log.add_node(prev_nd)
    print('done')
    return update_log.flush()