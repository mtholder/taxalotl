#!/usr/bin/env python
from __future__ import print_function

import io
import os
import shutil
import urllib
import urllib.request

from peyotl import (get_logger, write_as_json)

from taxalotl.ott_schema import (INP_OTT_SYNONYMS_HEADER,
                                 INP_OTT_TAXONOMY_HEADER,
                                 partition_ott_by_root_id)
from taxalotl.newick import normalize_newick
from taxalotl.partitions import (find_partition_dirs_for_taxonomy,
                                 has_any_partition_dirs,
                                 get_auto_gen_part_mapper,
                                 get_inp_taxdir,
                                 get_misc_inp_taxdir,
                                 get_taxon_partition, )
from taxalotl.tax_partition import TAX_SLICE_CACHE
from taxalotl.util import unlink, OutFile, OutDir

_LOG = get_logger(__name__)

def serialize_triple_object(o):
    return o.canonical_id if isinstance(o, SemGraphNode) else o

class SemGraphNode(object):
    def __init__(self, sem_graph, canoncial_id):
        self.canonical_id = canoncial_id
        self.graph = sem_graph
    def as_dict(self):
        d = {}
        for att in self.predicates:
            val = getattr(self, att, [])
            if val:
                d[att] = [serialize_triple_object(i) for i in val]
    @property
    def predicates(self):
        return []

def canonicalize(res_id, pred_id, entity_id):
    return '{}:{}:{}'.format(res_id, pred_id, entity_id)

class TaxonConceptSemNode(SemGraphNode):
    def __init__(self, sem_graph, canoncial_id):
        super(TaxonConceptSemNode, self).__init__(sem_graph, canoncial_id)
        self.is_child_of = None
        self.rank = None

    def claim_is_child_of(self, par_sem_node):
        self.is_child_of = par_sem_node

    def claim_rank(self, rank):
        self.rank = rank

    def predicates(self):
        return ['is_child_of', 'rank']

class SemGraph(object):
    att_list = ['_specimens', '_spec_based_names', '_group_names',
                '_taxon_concepts', '_authorities', '_references']
    att_set = frozenset(att_list)
    def __init__(self):
        for att in SemGraph.att_list:
            setattr(self, att, None)

    def add_taxon_concept(self, res_id, concept_id):
        return TaxonConceptSemNode(self, canonicalize(res_id, 'tc', concept_id))

    def __getattr__(self, item):
        hidden = '_{}'.format(item)
        if hidden not in SemGraph.att_set:
            raise AttributeError("'SemGraph' object has no attribute '{}'".format(item))
        v = getattr(self, hidden)
        if v is None:
            v = []
            setattr(self, hidden, v)
        return v

    def as_dict(self):
        d = {}
        for hidden in SemGraph.att_list:
            v = getattr(self, hidden)
            if v is not None:
                d[hidden[1:]] = {i.canonical_id(): i.as_dict() for i in v}
        return d


def semanticize_and_serialize_tax_part(taxolotl_config, res, fragment, out_dir, tax_part, tax_forest):
    sem_graph = semanticize_tax_part(taxolotl_config, res, fragment, tax_part, tax_forest)
    serialize_sem_graph(taxolotl_config, sem_graph, out_dir)

def semanticize_tax_part(taxolotl_config, res, fragment, tax_part, tax_forest):
    sem_graph = SemGraph()
    for r in tax_forest.roots.values():
        semanticize_subtree(sem_graph, res, r.root, par_sem_node=None)
    return sem_graph

def semanticize_subtree(sem_graph, res, node, par_sem_node=None):
    sem_node = res.semanticize_node_entry(sem_graph, node, par_sem_node)
    cr = node.children_refs
    csn = []
    if cr:
        csn = [semanticize_subtree(sem_graph, res, c, par_sem_node=sem_node) for c in cr]
    res.semanticize_node_exit(sem_graph, node, sem_node, child_sem_nodes=csn)
    return sem_node


def serialize_sem_graph(taxolotl_config, sem_graph, out_dir):
    with OutDir(out_dir):
        fp = os.path.join(out_dir, "sem_graph.json")
        with OutFile(fp) as out:
            write_as_json(sem_graph.as_dict(), out)