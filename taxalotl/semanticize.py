#!/usr/bin/env python
from __future__ import print_function

import os

from peyotl import (get_logger, write_as_json)

from taxalotl.util import OutFile, OutDir
from .taxonomic_ranks import (ABOVE_GENUS_SORTING_NUMBER,
                              SPECIES_SORTING_NUMBER,
                              GENUS_SORTING_NUMBER,
                              MINIMUM_HIGHER_TAXON_NUMBER)

from .name_parsing import (parse_genus_group_name,
                           parse_higher_name,
                           parse_name_string_without_context,
                           parse_sp_name, )
from .tax_partition import IGNORE_COMMON_NAME_SYN_TYPES
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
        self.problematic_synonyms = None
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
        return ['is_child_of', 'rank', 'has_name', 'id', 'undescribed', 'problematic_synonyms']

    @property
    def valid_combination(self):
        return None if self.has_name is None else self.has_name.valid_combination

    def claim_problematic_synonym_statement(self, name, syn_type, error_str):
        if self.problematic_synonyms is None:
            self.problematic_synonyms = []
        blob = {"name": name, "syn_type": syn_type, "problem": error_str}
        self.problematic_synonyms.append(blob)

    def claim_type_material(self, type_str):
        n = self.name_attached_to_type_specimen
        try:
            n.claim_type_material(type_str)
        except:
            e = 'could not find name that was type-material-based attached to TaxonConcept'
            self.claim_problematic_synonym_statement(type_str, 'type material', e)

    @property
    def name_attached_to_type_specimen(self):
        return None if self.has_name is None else self.has_name.name_attached_to_type_specimen

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
    extra_pred = ('combination', 'higher_group_name', 'genus_name',
                'subgenus_names', 'sp_epithet', 'infra_epithets', 'specimen_codes')
    name_sem_nd_pred = tuple(list(NameSemNode.name_sem_nd_pred) + list(extra_pred))

    def __init__(self, sem_graph, res_id, concept_id, name):
        super(VerbatimSemNode, self).__init__(sem_graph, res_id, 'combin', concept_id, name)
        self.combination = None
        self.higher_group_name = None
        self.genus_name = None
        self.subgenus_names = None
        self.sp_epithet = None
        self.infra_epithets = None
        self.specimen_codes = None

    @property
    def predicates(self):
        return VerbatimSemNode.name_sem_nd_pred

    def claim_higher_group_name(self, n):
        assert self.higher_group_name is None
        self.higher_group_name = n

    def claim_combination(self, n):
        assert self.combination is None
        self.combination = n

    def claim_genus(self, n):
        assert self.genus_name is None
        self.genus_name = n

    def claim_subgenus(self, n):
        if self.subgenus_names is None:
            self.subgenus_names = [n]
        else:
            self.subgenus_names.append(n)

    def claim_sp_epithet(self, n):
        assert self.sp_epithet is None
        self.sp_epithet = n

    def claim_infra_epithet(self, n):
        if self.infra_epithets is None:
            self.infra_epithets = [n]
        else:
            self.infra_epithets.append(n)

    def claim_specimen_code(self, n):
        if self.specimen_codes is None:
            self.specimen_codes = [n]
        else:
            self.specimen_codes.append(n)

    @property
    def valid_combination(self):
        return None if self.combination is None else self.combination.name

    @property
    def most_terminal_infra_epithet(self):
        if self.infra_epithets:
            if len(self.infra_epithets) > 1:
                x = [(len(i.name), i.name, i) for i in self.infra_epithets]
                x.sort()
                return x[-1][-1]
            return self.infra_epithets[0]
        return None

    @property
    def name_attached_to_type_specimen(self):
        ie = self.most_terminal_infra_epithet
        if ie is not None:
            return ie
        return self.sp_epithet

class GenusGroupSemNode(NameSemNode):
    def __init__(self, sem_graph, res_id, concept_id, name):
        super(GenusGroupSemNode, self).__init__(sem_graph, res_id, 'gen', concept_id, name)


class SpeciesGroupSemNode(NameSemNode):
    sp_grp_name_sem_nd_pred = tuple(list(NameSemNode.name_sem_nd_pred) + ['type_materials'])

    def __init__(self, sem_graph, res_id, concept_id, name):
        super(SpeciesGroupSemNode, self).__init__(sem_graph, res_id, 'sp', concept_id, name)
        self.type_materials = None

    def claim_type_material(self, type_str):
        if self.type_materials is None:
            self.type_materials = []
        self.type_materials.append(type_str)

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

def semanticize_node_synonym(res, sem_graph, node, sem_node, syn):
    try:
        rsn = node.rank_sorting_number()
    except KeyError:
        rsn = ABOVE_GENUS_SORTING_NUMBER
    if syn.syn_type == 'type material':
        if rsn >= MINIMUM_HIGHER_TAXON_NUMBER:
            _LOG.warn('PROBLEM: "{}" rsn={} node.rank={}'.format(syn.name, rsn, node.rank))
            sem_node.claim_problematic_synonym_statement(syn.name, syn.syn_type, "type_material for higher taxon")
        else:
            sem_node.claim_type_material(syn.name.strip())
        return
    name_dict = parse_name_string_without_context(syn.name)
    if not name_dict:
        _LOG.debug('WARNING: Could not parse "{}" as a {} for {} ({})'.format(syn.name, syn.syn_type, node.name, node.id))
        sem_node.claim_problematic_synonym_statement(syn.name, syn.syn_type, "unparseable")
    '''

    

    else:
        valid_combo = sem_node.valid_combination
        if len(syn.name.split()) != len(valid_combo.split()):
            _LOG.debug('WARNING: "{}" is a {} for {} ({})'.format(syn.name, syn.syn_type, node.name, node.id))
        else:
            _LOG.debug('         "{}" is a {} for {} ({})'.format(syn.name, syn.syn_type, node.name, node.id))
    st = syn.syn_type
    '''

def semanticize_node_auth_synonym(res, sem_graph, node, sem_node, syn):
    _LOG.debug('"{}" is a {} for {} ({})'.format(syn.name, syn.syn_type, node.name, node.id))
    st = syn.syn_type

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
    # First, Preorder sweep of the tree to deal with taxon concepts and valid names
    nd_sem_pairs = []
    for r in tax_forest.roots.values():
        semanticize_subtree(sem_graph, res, r.root, nd_sem_pairs, par_sem_node=None)
    # grab all the synonyms
    t2syn = tax_part.parsed_synonyms_by_id(ignored_syn_types=IGNORE_COMMON_NAME_SYN_TYPES)
    delay_auth_syns = []
    # process the non-"authority" synonyms next
    for node, sem_node in nd_sem_pairs:
        if sem_node is not None:
            for syn in t2syn.get(node.id, []):
                if syn.syn_type == 'authority':
                    delay_auth_syns.append((node, sem_node, syn))
                    continue
                res.semanticize_node_synonyms(sem_graph, node, sem_node, syn)
    # now the "authority" synonyms...
    for node, sem_node, syn in delay_auth_syns:
        res.semanticize_node_authority_synonyms(sem_graph, node, sem_node, syn)
    return sem_graph


def semanticize_subtree(sem_graph, res, node, nd_sem_pairs, par_sem_node=None):
    sem_node = res.semanticize_node_entry(sem_graph, node, par_sem_node)
    nd_sem_pairs.append((node, sem_node))
    if sem_node is None:
        return None
    cr = node.children_refs
    csn = []
    if cr:
        csn = [semanticize_subtree(sem_graph, res, c,
                                   nd_sem_pairs, par_sem_node=sem_node) for c in cr]
    csn = [i for i in csn if i is not None]
    res.semanticize_node_exit(sem_graph, node, sem_node, child_sem_nodes=csn)
    return sem_node


# noinspection PyUnusedLocal
def serialize_sem_graph(taxolotl_config, sem_graph, out_dir):
    with OutDir(out_dir):
        fp = os.path.join(out_dir, "sem_graph.json")
        with OutFile(fp) as out:
            write_as_json(sem_graph.as_dict(), out, indent=2)
