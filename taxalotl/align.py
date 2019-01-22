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


def _get_findable_names(leaf, tree):
    r = []
    if leaf.best_rank_sort_number > SPECIES_SORTING_NUMBER:
        return r
    findability_score = 10*len(leaf.src_dict) - len(leaf.synonyms)
    r.append((findability_score, id(leaf), leaf))
    while leaf.best_rank_sort_number < SPECIES_SORTING_NUMBER:
        try:
            leaf = tree.id_to_taxon[leaf.par_id]
        except:
            break
        if leaf.best_rank_sort_number <= SPECIES_SORTING_NUMBER:
            findability_score = 10 * len(leaf.src_dict) - len(leaf.synonyms)
            r.append((findability_score, id(leaf), leaf))
    return r


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
    scored_leaves = []
    for n, leaf in enumerate(ott_tree.leaves()):
        scored_leaves.extend(_get_findable_names(leaf, ott_tree))
    scored_leaves.sort(reverse=True)
    ott_leaf_label_list = []
    ott_lls = set()
    for i in scored_leaves:
        n = i[2].name
        if n not in ott_lls:
            ott_leaf_label_list.append(n)
            ott_lls.add(n)
    _LOG.info('Will look for {} <= species taxa names...'.format(len(ott_leaf_label_list)))
    res_forest = res.get_taxon_forest_for_partition(part_name)
    if res_forest:
        _LOG.info('{} already separated for {}'.format(res.id, part_name))
        assert len(res_forest.trees) == 1
        rt = res_forest.trees[0]
    else:
        rt = separate_based_on_tip_overlap(taxalotl_config, ott_res, ott_lls, ott_tree, res, part_name)
    sys.exit('rt = {}'.format(rt))

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
        res_syn = res.get_parsed_synonyms_by_id(higher_part_name)
        slice_tree.attach_parsed_synonyms_set(res_syn)
        slice_tree.add_best_guess_rank_sort_number()
        n = 0
        for n, nd in enumerate(slice_tree.postorder()):
            if nd.children_refs:
                nd.found_names, nd.unfound_names = set(), set()
                for c in nd.children_refs:
                    nd.found_names.update(c.found_names)
                    nd.unfound_names.update(c.unfound_names)
            else:
                nd.found_names = set()
                nd.unfound_names = set()
                if nd.name in ott_lls:
                    nd.found_names.add(nd.name)
                    # _LOG.debug('found "{}"'.format(nd.name))
                else:
                    nd.unfound_names.add(nd.name)
        root_found = slice_tree.root.found_names
        pf = tot_leaves.intersection(root_found)
        if pf:
            m = 'Leaves {} found in multiple trees for {} at {}'
            raise ValueError(m.format(pf, res.id, higher_part_name))
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



