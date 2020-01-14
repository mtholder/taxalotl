#!/usr/bin/env python
from __future__ import print_function

import os

from peyotl import (get_logger, write_as_json)

from taxalotl.util import OutFile, OutDir
from .taxonomic_ranks import SPECIES_SORTING_NUMBER, GENUS_SORTING_NUMBER
from .name_parsing import (parse_genus_group_name, parse_higher_name, parse_sp_name)

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
                if isinstance(val, list) or isinstance(val, tuple):
                    d[att] = [serialize_triple_object(i) for i in val]
                else:
                    d[att] = serialize_triple_object(val)
        return d

    @property
    def predicates(self):
        return []


def canonicalize(res_id, pred_id, entity_id):
    return '{}:{}:{}'.format(res_id, pred_id, entity_id)


class TaxonConceptSemNode(SemGraphNode):
    def __init__(self, sem_graph, res, concept_id):
        res_id = res.base_resource.id
        self.concept_id = concept_id
        ci = canonicalize(res_id, 'tc', concept_id)
        super(TaxonConceptSemNode, self).__init__(sem_graph, ci)
        self.is_child_of = None
        self.rank = None
        self.has_name = None
        self.undescribed = None
        self.id = '{}:{}'.format(res_id, concept_id)

    def claim_is_child_of(self, par_sem_node):
        assert self.is_child_of is None
        self.is_child_of = par_sem_node

    def claim_rank(self, rank):
        assert self.rank is None
        self.rank = rank

    def claim_name(self, name_sem_node):
        assert self.has_name is None
        self.has_name = name_sem_node

    def claim_undescribed(self):
        self.undescribed = True

    @property
    def predicates(self):
        return ['is_child_of', 'rank', 'has_name', 'id', 'undescribed']


class NameSemNode(SemGraphNode):
    name_sem_nd_pred = ('name', )
    def __init__(self, sem_graph, res_id, tag, concept_id, name):
        ci = canonicalize(res_id, tag, concept_id)
        super(NameSemNode, self).__init__(sem_graph, ci)
        self._name = name

    @property
    def name(self):
        return self._name

    @property
    def predicates(self):
        return NameSemNode.name_sem_nd_pred


class CombinationSemNode(NameSemNode):
    def __init__(self, sem_graph, res_id, concept_id, name):
        super(CombinationSemNode, self).__init__(sem_graph, res_id, 'combin', concept_id, name)


class VerbatimSemNode(NameSemNode):
    extra_pred = ('combinations', 'higher_group_names', 'genus_names',
                'subgenus_names', 'sp_epithets', 'infra_epithets')
    name_sem_nd_pred = tuple(list(NameSemNode.name_sem_nd_pred) + list(extra_pred))

    def __init__(self, sem_graph, res_id, concept_id, name):
        super(VerbatimSemNode, self).__init__(sem_graph, res_id, 'combin', concept_id, name)
        self.combinations = []
        self.higher_group_names = []
        self.genus_names = []
        self.subgenus_names = []
        self.sp_epithets = []
        self.infra_epithets = []
        self.specimen_codes = []

    @property
    def predicates(self):
        return VerbatimSemNode.name_sem_nd_pred

    def claim_higher_group_name(self, n):
        self.higher_group_names.append(n)

    def claim_combination(self, n):
        self.combinations.append(n)

    def claim_genus(self, n):
        self.genus_names.append(n)

    def claim_subgenus(self, n):
        self.subgenus_names.append(n)

    def claim_sp_epithet(self, n):
        self.sp_epithets.append(n)

    def claim_infra_epithet(self, n):
        self.infra_epithets.append(n)

    def claim_specimen_code(self, n):
        self.specimen_codes.append(n)


class GenusGroupSemNode(NameSemNode):
    def __init__(self, sem_graph, res_id, concept_id, name):
        super(GenusGroupSemNode, self).__init__(sem_graph, res_id, 'gen', concept_id, name)


class SpeciesGroupSemNode(NameSemNode):
    def __init__(self, sem_graph, res_id, concept_id, name):
        super(SpeciesGroupSemNode, self).__init__(sem_graph, res_id, 'sp', concept_id, name)


class HigherGroupSemNode(NameSemNode):
    def __init__(self, sem_graph, res_id, concept_id, name):
        super(HigherGroupSemNode, self).__init__(sem_graph, res_id, 'clade', concept_id, name)


class SpecimenCodeSemNode(NameSemNode):
    def __init__(self, sem_graph, res_id, concept_id, name):
        super(SpecimenCodeSemNode, self).__init__(sem_graph, res_id, 'spec_code', concept_id, name)


def _find_by_name(container, name):
    if container:
        for el in container:
            if el._name == name:
                return el
    return None


class SemGraph(object):
    att_list = ['_specimens',
                '_specimen_codes',
                '_species_group_epithets',
                '_genus_group_names',
                '_combinations',
                '_higher_group_names',
                '_verbatim_name',
                '_taxon_concepts', '_authorities', '_references']
    att_set = frozenset(att_list)

    def __init__(self):
        for att in SemGraph.att_list:
            setattr(self, att, None)

    def _add_name(self, container, node_type, res_id, concept_id, name):
        x = _find_by_name(container, name)
        if x is None:
            x = node_type(self, res_id, concept_id, name)
            container.append(x)
        return x

    def add_combination(self, res_id, concept_id, name):
        return self._add_name(self.combinations, CombinationSemNode, res_id, concept_id, name)

    def add_verbatim_name(self, res_id, concept_id, name):
        return self._add_name(self.verbatim_name, VerbatimSemNode, res_id, concept_id, name)

    def add_genus(self, res_id, concept_id, name):
        return self._add_name(self.genus_group_names, GenusGroupSemNode, res_id, concept_id, name)

    add_subgenus = add_genus

    def add_higher_group_name(self, res_id, concept_id, name):
        return self._add_name(self.higher_group_names, HigherGroupSemNode,
                              res_id, concept_id, name)

    def add_sp_epithet(self, res_id, concept_id, name):
        return self._add_name(self.species_group_epithets, SpeciesGroupSemNode,
                              res_id, concept_id, name)

    add_infra_epithet = add_sp_epithet

    def add_specimen_code(self, res_id, concept_id, name):
        return self._add_name(self.specimen_codes, SpecimenCodeSemNode,
                              res_id, concept_id, name)

    def add_taxon_concept(self, res, concept_id):
        x = TaxonConceptSemNode(self, res, concept_id)
        self.taxon_concepts.append(x)
        return x

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
                d[hidden[1:]] = {i.canonical_id: i.as_dict() for i in v}
        return d


def semanticize_node_name(res, sem_graph, tc, node):
    name_dict = None
    try:
        rsn = node.rank_sorting_number()
    except KeyError:
        pass
    else:
        if rsn is None:
            pass
        elif rsn <= SPECIES_SORTING_NUMBER:
            name_dict = parse_sp_name(node.name, node.rank)
        elif rsn <= GENUS_SORTING_NUMBER:
            name_dict = parse_genus_group_name(node.name, node.rank)
    if name_dict is None:
        name_dict = parse_higher_name(node.name, node.rank)
    semanticize_names(res, sem_graph, tc, node.name, name_dict)


def semanticize_names(res, sem_graph, taxon_concept_sem_node, name, name_dict):
    tcsn = taxon_concept_sem_node
    if name_dict.get('undescribed'):
        tcsn.claim_undescribed()
    rn = sem_graph.add_verbatim_name(res.base_resource.id, tcsn.concept_id, name)
    name_part_holder = rn
    tcsn.claim_name(rn)
    combination = name_dict.get('combination')
    bresid = res.base_resource.id
    if combination:
        cn = sem_graph.add_combination(bresid, tcsn.concept_id, combination)
        rn.claim_combination(cn)
    genus = name_dict.get('genus')
    if genus:
        cn = sem_graph.add_genus(bresid, tcsn.concept_id, genus)
        name_part_holder.claim_genus(cn)
    subgenus = name_dict.get('subgenus')
    if subgenus:
        cn = sem_graph.add_subgenus(bresid, tcsn.concept_id, subgenus)
        name_part_holder.claim_subgenus(cn)
    sp_epithet = name_dict.get('sp_epithet')
    if sp_epithet:
        cn = sem_graph.add_sp_epithet(bresid, tcsn.concept_id, sp_epithet)
        name_part_holder.claim_sp_epithet(cn)
    infra_epithet = name_dict.get('infra_epithet')
    if infra_epithet:
        cn = sem_graph.add_infra_epithet(bresid, tcsn.concept_id, infra_epithet)
        name_part_holder.claim_infra_epithet(cn)
    higher_group_name = name_dict.get('higher_group_name')
    if higher_group_name:
        cn = sem_graph.add_higher_group_name(bresid, tcsn.concept_id, higher_group_name)
        name_part_holder.claim_higher_group_name(cn)
    specimen_code = name_dict.get('specimen_code')
    if specimen_code:
        cn = sem_graph.add_specimen_code(bresid, tcsn.concept_id, specimen_code)
        name_part_holder.claim_specimen_code(cn)
    return rn


def semanticize_and_serialize_tax_part(taxolotl_config, res, fragment, out_dir, tax_part, tax_forest):
    sem_graph = semanticize_tax_part(taxolotl_config, res, fragment, tax_part, tax_forest)
    serialize_sem_graph(taxolotl_config, sem_graph, out_dir)


# noinspection PyUnusedLocal
def semanticize_tax_part(taxolotl_config, res, fragment, tax_part, tax_forest):
    sem_graph = SemGraph()
    for r in tax_forest.roots.values():
        semanticize_subtree(sem_graph, res, r.root, par_sem_node=None)
    return sem_graph


def semanticize_subtree(sem_graph, res, node, par_sem_node=None):
    sem_node = res.semanticize_node_entry(sem_graph, node, par_sem_node)
    if sem_node is None:
        return None
    cr = node.children_refs
    csn = []
    if cr:
        csn = [semanticize_subtree(sem_graph, res, c, par_sem_node=sem_node) for c in cr]
    csn = [i for i in csn if i is not None]
    res.semanticize_node_exit(sem_graph, node, sem_node, child_sem_nodes=csn)
    return sem_node


# noinspection PyUnusedLocal
def serialize_sem_graph(taxolotl_config, sem_graph, out_dir):
    with OutDir(out_dir):
        fp = os.path.join(out_dir, "sem_graph.json")
        with OutFile(fp) as out:
            write_as_json(sem_graph.as_dict(), out, indent=2)
