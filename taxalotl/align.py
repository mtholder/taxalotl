#!/usr/bin/env python
from __future__ import print_function

import sys
import os
from enum import IntEnum, IntFlag
from typing import Dict, List
import copy

from peyotl import (get_logger, write_as_json)

from .config import TaxalotlConfig
from .tree import TaxonTree
from .partitions import (PART_NAMES)
from .resource_wrapper import TaxonomyWrapper
from .taxonomic_ranks import (GENUS_RANK_TO_SORTING_NUMBER,
                              MAX_INFRASPECIFIC_NUMBER,
                              MINIMUM_HIGHER_TAXON_NUMBER,
                              SPECIES_SORTING_NUMBER)
from .tax_partition import get_taxon_partition
from .util import get_true_false_repsonse

_LOG = get_logger(__name__)
out_stream = sys.stdout

class NodeFilter(IntEnum):
    SPECIES = 1
    SP_OR_BELOW = 2
    TIP = 3

def align_resource(taxalotl_config: TaxalotlConfig,
                   ott_res: TaxonomyWrapper,
                   res: TaxonomyWrapper,
                   level_list: List[str]):
    m = 'Could not align taxonomy because {} has not been partitioned.'
    for el in [ott_res, res]:
        if not el.has_been_partitioned():
            raise RuntimeError(m.format(el.id))
    if level_list == [None]:
        level_list = PART_NAMES
    for part_name in level_list:
        align_for_level(taxalotl_config, ott_res, res, part_name)

def _register_name(tup_list, name_to_ott_id_list, leaf, name, ott_id):
    find_score = 10 * len(leaf.src_dict) - len(leaf.synonyms)
    tup_list.append((find_score, id(leaf), name))
    name_to_ott_id_list.setdefault(name, set()).add(ott_id)

def _get_findable_names(leaf, tree, name_to_ott_id_list, nd_filter, include_synonyms):
    if nd_filter == NodeFilter.SPECIES and leaf.best_rank_sort_number != SPECIES_SORTING_NUMBER:
        return []
    if leaf.best_rank_sort_number > SPECIES_SORTING_NUMBER and nd_filter == NodeFilter.SP_OR_BELOW:
        return []
    if nd_filter == NodeFilter.TIP and leaf.children_refs:
        return []
    r = []
    _register_name(r, name_to_ott_id_list, leaf, leaf.name, leaf.id)
    if include_synonyms:
        for syn in leaf.synonyms:
            _register_name(r, name_to_ott_id_list, leaf, syn.name, leaf.id)
    if nd_filter == NodeFilter.SP_OR_BELOW:
        while leaf.best_rank_sort_number < SPECIES_SORTING_NUMBER:
            try:
                leaf = tree.id_to_taxon[leaf.par_id]
            except:
                break
            if leaf.best_rank_sort_number <= SPECIES_SORTING_NUMBER:
                _register_name(r, name_to_ott_id_list, leaf, leaf.name, leaf.id)
                if include_synonyms:
                    for syn in leaf.synonyms:
                        _register_name(r, name_to_ott_id_list, leaf, syn.name, leaf.id)
    return r

def get_findable_names(ott_tree, node_filter=NodeFilter.SP_OR_BELOW, include_synonyms=False):
    scored_nodes = []
    name_to_ott_id_set = {}
    gen = ott_tree.leaves if node_filter != NodeFilter.SPECIES else ott_tree.postorder
    for node in gen():
        x = _get_findable_names(node, ott_tree, name_to_ott_id_set, node_filter, include_synonyms)
        scored_nodes.extend(x)
    scored_nodes.sort(reverse=True)
    ott_leaf_label_list = []
    ott_lls = set()
    for i in scored_nodes:
        n = i[2]
        if n not in ott_lls:
            ott_leaf_label_list.append(n)
            ott_lls.add(n)
    return ott_leaf_label_list, ott_lls, name_to_ott_id_set

def align_for_level(taxalotl_config: TaxalotlConfig,
                    ott_res: TaxonomyWrapper,
                    res: TaxonomyWrapper,
                    part_name: str):
    fragment = taxalotl_config.get_fragment_from_part_name(part_name)
    _LOG.info('align for {} for {}'.format(fragment, res.id))
    ott_forest = ott_res.get_taxon_forest_for_partition(part_name)
    assert len(ott_forest.trees) == 1
    ott_tree = ott_forest.trees[0]
    prev_syn = ott_res.get_parsed_synonyms_by_id(part_name)
    ott_tree.attach_parsed_synonyms_set(prev_syn)
    ott_tree.add_best_guess_rank_sort_number()

    ott_leaf_label_list, ott_lls, name_to_ott_id_set = get_findable_names(ott_tree)

    _LOG.info('Will look for {} <= species taxa names...'.format(len(ott_leaf_label_list)))
    res_forest = res.get_taxon_forest_for_partition(part_name)
    if res_forest:
        _LOG.info('{} already separated for {}'.format(res.id, part_name))
    else:
        separate_based_on_tip_overlap(taxalotl_config, ott_res, ott_lls, ott_tree, res, part_name)
        res_forest = res.get_taxon_forest_for_partition(part_name)
        if not res_forest:
            m = 'Failed to separate {} for {}'
            raise ValueError(m.format(res.id, part_name))
    non_incert_trees, incert_trees = [], []
    for tree in res_forest.trees:
        if tree.root.flags and 'incertae_sedis' in tree.root.flags:
            incert_trees.append(tree)
        else:
            non_incert_trees.append(tree)
    align_trees_for_level(ott_res, ott_tree, res, part_name, non_incert_trees, incert_trees)

def align_trees_for_level(ott_res, ott_tree, res, part_name, non_incert_trees, incert_trees):
    assert len(non_incert_trees) == 1 # should be NotImplementedError
    sp_tup = get_findable_names(ott_tree, node_filter=NodeFilter.SPECIES, include_synonyms=False)
    sp_ott_ls = sp_tup[1]
    print('non_incert_trees =', non_incert_trees)
    print('incert_trees =', incert_trees)
    tree_l = non_incert_trees + incert_trees
    for tree in tree_l:
        attach_synonyms_and_find_strict_name_matches(res, tree, part_name, sp_ott_ls)
        for nd in tree.postorder():
            nd.match_status = None
    found_tip_names = {}

    _new_match_stat(tree_l, sp_ott_ls, MatchStatus.VALID_SP_OTT_VALID_SP, NodeFilter.SPECIES, False)
    _new_match_stat(tree_l, sp_ott_ls, MatchStatus.SYN_SP_OTT_VALID_SP, NodeFilter.SPECIES, True)
    _new_match_stat(tree_l, sp_ott_ls, MatchStatus.VALID_INF_OTT_VALID_SP, NodeFilter.SP_OR_BELOW, False)
    _new_match_stat(tree_l, sp_ott_ls, MatchStatus.SYN_INF_OTT_VALID_SP, NodeFilter.SP_OR_BELOW, True)
    infsp_tup = get_findable_names(ott_tree, node_filter=NodeFilter.SP_OR_BELOW, include_synonyms=False)
    infsp_ott_ls = infsp_tup[1]
    _new_match_stat(tree_l, infsp_ott_ls, MatchStatus.VALID_SP_OTT_VALID_INF, NodeFilter.SPECIES, False)
    _new_match_stat(tree_l, infsp_ott_ls, MatchStatus.SYN_SP_OTT_VALID_INF, NodeFilter.SPECIES, True)
    _new_match_stat(tree_l, infsp_ott_ls, MatchStatus.VALID_INF_OTT_VALID_INF, NodeFilter.SP_OR_BELOW,
                    False)
    _new_match_stat(tree_l, infsp_ott_ls, MatchStatus.SYN_INF_OTT_VALID_INF, NodeFilter.SP_OR_BELOW, True)
    for tree in tree_l:
        for nd in tree.postorder():
            if nd.best_rank_sort_number <= SPECIES_SORTING_NUMBER and nd.match_status is None:
                print('Still unmatched: {}'.format(str(nd)))

def _new_match_stat(tree_l, ott_lls, match_stat, nd_filter, check_synonyms):
    for tree in tree_l:
        mark_found_unfound_name_matches(tree,
                                        ott_lls,
                                        nd_filter=nd_filter,
                                        check_synonyms=check_synonyms)
        for nd in tree.postorder():
            if nd.match_status is None and nd.matched_to_name is not None:
                nd.match_status = match_stat
                if nd.name != nd.matched_to_name:
                    m = '{} match for "{}" (valid = "{}")'
                    print(m.format(match_stat.name, nd.matched_to_name, nd.name))
                else:
                    m = '{} match for "{}"'
                    print(m.format(match_stat.name, nd.matched_to_name))


class MatchStatus(IntFlag):
    EXT_VALID =  0x001
    EXT_SYN =    0x002
    OTT_VALID =  0x004
    OTT_SYN =    0x008
    EXT_SP =     0x010
    EXT_INF_SP = 0x020
    EXT_CLADE =  0x040
    OTT_SP =     0x080
    OTT_INF_SP = 0x100
    OTT_CLADE =  0x200
    VALID_SP_OTT_VALID_SP = EXT_VALID | OTT_VALID | EXT_SP | OTT_SP
    SYN_SP_OTT_VALID_SP = EXT_SYN | OTT_VALID | EXT_SP | OTT_SP
    VALID_INF_OTT_VALID_SP = EXT_VALID | OTT_VALID | EXT_INF_SP | OTT_SP
    SYN_INF_OTT_VALID_SP = EXT_SYN | OTT_VALID | EXT_INF_SP | OTT_SP
    VALID_SP_OTT_VALID_INF = EXT_VALID | OTT_VALID | EXT_SP | OTT_INF_SP
    SYN_SP_OTT_VALID_INF = EXT_SYN | OTT_VALID | EXT_SP | OTT_INF_SP
    VALID_INF_OTT_VALID_INF = EXT_VALID | OTT_VALID | EXT_INF_SP | OTT_INF_SP
    SYN_INF_OTT_VALID_INF = EXT_SYN | OTT_VALID | EXT_INF_SP | OTT_INF_SP


def attach_synonyms_and_find_strict_name_matches(res, tree, part_name, ott_lls):
    res_syn = res.get_parsed_synonyms_by_id(part_name)
    tree.attach_parsed_synonyms_set(res_syn, warn_missing_target=False)
    tree.add_best_guess_rank_sort_number()
    mark_found_unfound_name_matches(tree, ott_lls)


def mark_found_unfound_name_matches(tree,
                                    ott_lls,
                                    nd_filter=NodeFilter.TIP,
                                    check_synonyms=False):
    for nd in tree.postorder():
        if nd_filter == NodeFilter.TIP:
            do_name_check = not nd.children_refs
            do_des_union = not do_name_check
        elif nd_filter == NodeFilter.SPECIES:
            do_name_check = nd.best_rank_sort_number == SPECIES_SORTING_NUMBER
            do_des_union = nd.best_rank_sort_number > SPECIES_SORTING_NUMBER
        else:
            assert nd_filter == NodeFilter.SP_OR_BELOW
            do_name_check = nd.best_rank_sort_number <= SPECIES_SORTING_NUMBER
            do_des_union = bool(nd.children_refs)
        nd.found_names, nd.unfound_names = set(), set()
        nd.matched_to_name = None
        if do_name_check:
            found = False
            if nd.name in ott_lls:
                nd.found_names.add(nd.name)
                nd.matched_to_name = nd.name
            else:
                if check_synonyms:
                    for syn in nd.synonyms:
                        if syn.name in ott_lls:
                            nd.found_names.add(syn.name)
                            nd.matched_to_name = syn.name
                            break
            if nd.matched_to_name is None:
                nd.unfound_names.add(nd.name)
        if do_des_union:
            if nd.children_refs:
                for c in nd.children_refs:
                    nd.found_names.update(c.found_names)
                    nd.unfound_names.update(c.unfound_names)


def separate_based_on_tip_overlap(taxalotl_config, ott_res, ott_lls, ott_tree, res, part_name):
    fragment = taxalotl_config.get_fragment_from_part_name(part_name)
    orig_fragment = fragment
    higher_part_name = part_name
    res_forest = None
    while not res_forest:
        fragment = os.path.split(fragment)[0]
        higher_part_name = os.path.split(fragment)[-1]
        if higher_part_name == 'partitioned':
            raise ValueError('Could not find any parition for {}'.format(res.id))
        res_forest = res.get_taxon_forest_for_partition(higher_part_name)
        if res_forest:
            _LOG.info('{} trees for {} at {}'.format(len(res_forest.trees), res.id, higher_part_name))
        else:
            _LOG.info('no trees for {} at {}'.format(res.id, higher_part_name))
    tot_leaves = set()
    to_move_ids = set()
    for tree_ind, slice_tree in enumerate(res_forest.trees):
        attach_synonyms_and_find_strict_name_matches(res, slice_tree, higher_part_name, ott_lls)
        root_found = slice_tree.root.found_names
        pf = tot_leaves.intersection(root_found)
        if pf:
            m = 'Leaves {} found in multiple trees for {} at {}'
            raise ValueError(m.format(pf, res.id, part_name))
        tot_leaves.update(root_found)
        m = 'Found {}/{} names in tree_ind={} for {} at {}'
        _LOG.info(m.format(len(root_found), len(ott_lls), tree_ind, res.id, higher_part_name))

        only_overlap = []
        mainly_overlap = []
        if len(root_found) > 0:
            curr_node = slice_tree.root
            curr_found_names = copy.copy(curr_node.found_names)
            while True:
                next_node = None
                some_overlap = []
                for c in curr_node.children_refs:
                    if c.found_names == curr_found_names:
                        _LOG.info('Moving tipward from {} to {}'.format(curr_node.name, c.name))
                        next_node = c
                        break
                    else:
                        if c.found_names:
                            some_overlap.append(c)
                            if not c.unfound_names:
                                only_overlap.append(c)
                            elif len(c.found_names) > len(c.unfound_names):
                                mainly_overlap.append(c)
                            m = '  {} has {} relevant names and {} irrelevant'
                            _LOG.info(m.format(c.name, len(c.found_names), len(c.unfound_names)))
                if next_node:
                    curr_node = next_node
                else:
                    break
            _LOG.info('MRCA of tips is {}'.format(curr_node.name))
            to_move = only_overlap + mainly_overlap
            if not to_move:
                to_move_ids.add(curr_node.id)
            else:
                mtmplate = '  moving {} ({} relevant/ {} irrelevant) names;'
                msg_list = ['Perform separation by']
                for node in to_move:
                    m = mtmplate.format(node.name, len(node.found_names), len(node.unfound_names))
                    msg_list.append(m)
                msg_list.append('?  Enter y to confirm:')
                prompt = '\n'.join(msg_list)
                if not get_true_false_repsonse(prompt, def_value=True):
                    return None
                to_move_ids.update({i.id for i in to_move})
    if not to_move_ids:
        return None
    from taxalotl.dynamic_partitioning import perform_dynamic_separation
    root_ott_taxon = ott_tree.root
    new_sep_val = {
        "name": root_ott_taxon.name,
        "uniqname": root_ott_taxon.name_that_is_unique,
        "src_dict": {res.base_id: list(to_move_ids)}
    }
    sbo = {root_ott_taxon.id: new_sep_val}
    perform_dynamic_separation(ott_res, res, higher_part_name, separation_by_ott=sbo)



