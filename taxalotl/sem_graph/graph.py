#!/usr/bin/env python
from __future__ import print_function

from peyotl import (get_logger, )
from .taxon_concept_node import TaxonConceptSemNode
from .verbatim_name import VerbatimSemNode
from .names_for_ranks import (SpeciesGroupSemNode,
                              GenusGroupSemNode,
                              HigherGroupSemNode,
                              SpecimenCodeSemNode,
                              )
from .name import CombinationSemNode
from .graph_node import AuthoritySemNode

_LOG = get_logger(__name__)


class SemGraph(object):
    att_list = ['_authorities',
                '_specimens',
                '_specimen_codes',
                '_species_group_epithets',
                '_genus_group_names',
                '_combinations',
                '_higher_group_names',
                '_verbatim_name',
                '_taxon_concepts', '_references',
                ]
    att_set = frozenset(att_list)

    def register_obj(self, id_minting_context, obj):
        return _register_new_node(self, id_minting_context, obj)

    def __init__(self, taxolotl_config, res):
        self.config = taxolotl_config
        self.res = res
        self._by_id = {}
        self._authorities = None
        self._specimens = None
        self._specimen_codes = None
        self._species_group_epithets = None
        self._genus_group_names = None
        self._combinations = None
        self._higher_group_names = None
        self._verbatim_name = None
        self._taxon_concepts = None
        self._references = None

    def denormalize_homonyms(self):
        # for species ranK:
        #   multiple authority entries
        #   same valid epithet in multiple valid genera
        to_mint = []
        for tc in self.taxon_concept_list:
            if tc.rank and tc.rank == 'species':
                epithet = tc.most_terminal_name
                if epithet is None:
                    _LOG.warn('NO Epithet for  = {}'.format(tc.__dict__))
                    continue
                if isinstance(epithet._authority, list):
                    _LOG.warn('epithet = {}'.format(epithet.__dict__))
        # import sys
        # sys.exit(1)

    @property
    def taxon_concept_list(self):
        return self._taxon_concepts if self._taxon_concepts else []

    def _all_specimen_based_tax_con_dict(self):
        r = {}
        for tcobj in self.taxon_concept_list:
            if tcobj.is_specimen_based:
                r[tcobj.canonical_id] = tcobj
        return r

    def _all_higher_tax_con_dict(self):
        r = {}
        for tcobj in self.taxon_concept_list:
            if not tcobj.is_specimen_based:
                r[tcobj.canonical_id] = tcobj
        return r

    def find_valid_genus(self, genus_name):
        r = []
        for tc in self._all_higher_tax_con_dict().values():
            if tc.rank and tc.rank == 'genus':
                if tc.is_valid_for_name(genus_name):
                    r.append(tc)
        return r

    def specimen_based_synonym_taxa(self):
        d = self._all_specimen_based_tax_con_dict()
        r = {}
        for tcid, tcobj in d.items():
            syn_list = tcobj.synonyms if tcobj.synonyms else []
            for syn_id in syn_list:
                r[syn_id] = d[syn_id]
        return r

    @property
    def valid_specimen_based_taxa(self):
        d = self._all_specimen_based_tax_con_dict()
        r = {}
        for tcid, tcobj in d.items():
            if not tcobj.is_synonym_of:
                r[tcid] = tcobj
        return r

    @property
    def valid_taxa_dict(self):
        raw = self._taxon_concepts if self._taxon_concepts else []
        r = {}
        for tcobj in raw:
            if not tcobj.is_synonym_of:
                r[tcobj.canonical_id] = tcobj
        return r

    @property
    def valid_name_to_taxon_concept_map(self):
        return {i.valid_name.name: i for i in self.valid_taxa_dict.values()}

    def get_by_id(self, can_id, default=None):
        return self._by_id.get(can_id, default)

    def _add_name(self, container, node_type, parent_sem_node, name, extra_container=None):
        search_cont = container if extra_container is None else extra_container.contained
        x = _find_by_name(search_cont, name)
        if x is None:
            d = {'parent_id': parent_sem_node.canonical_id}
            if extra_container is not None:
                d['class_tag'] = 'epi'  # need to figure out if this is the best choice for any extra container obj
            x = node_type(self, d, name)
            search_cont.append(x)
            if search_cont is not container:
                container.append(x)
        return x

    def add_authority(self, res_id, name_sem, authors, year):
        auth_list = self.authorities
        x = None
        for a in auth_list:
            if a.authors == authors and a.year == year:
                x = a
                break
        if x is None:
            d = {'parent_id': name_sem.canonical_id}
            x = AuthoritySemNode(self, d, authors, year)
        auth_list.append(x)
        name_sem.claim_authority(x)
        return x

    def _add_normalized(self, par_sem_node, name_str):
        return self._add_name(self.combinations, CombinationSemNode, par_sem_node, name_str)

    def _add_combination(self, par_sem_node, name_str):
        return self._add_name(self.combinations, CombinationSemNode, par_sem_node, name_str)

    def _add_verbatim_name(self, tax_con_sem_node, name_str):
        return self._add_name(self.verbatim_name, VerbatimSemNode, tax_con_sem_node, name_str)

    def _add_genus(self, par_sem_node, name_str):
        return self._add_name(self.genus_group_names, GenusGroupSemNode, par_sem_node, name_str)

    _add_subgenus = _add_genus

    def _add_higher_group_name(self, par_sem_node, name_str):
        return self._add_name(self.higher_group_names, HigherGroupSemNode, par_sem_node, name_str)

    def _add_sp_epithet(self, par_sem_node, name_str, prev_word_sn):
        return self._add_name(self.species_group_epithets, SpeciesGroupSemNode, par_sem_node, name_str, prev_word_sn)

    _add_infra_epithet = _add_sp_epithet

    def _add_specimen_code(self, par_sem_node, name_str):
        return self._add_name(self.specimen_codes, SpecimenCodeSemNode, par_sem_node, name_str)

    def add_taxon_concept(self, foreign_id):
        x = TaxonConceptSemNode(self, foreign_id)
        self.taxon_concepts.append(x)
        return x

    def __getattr__(self, item):
        hidden = '_{}'.format(item)
        if hidden not in SemGraph.att_set:
            return self.__getattribute__(item)
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


def _find_by_name(container, name):
    if container:
        for el in container:
            if el._name == name:
                return el
    return None


_NODE_CLASS_NAME_TO_CT = {AuthoritySemNode: 'auth',
                          SpecimenCodeSemNode: 'spec_code',
                          HigherGroupSemNode: 'clade',
                          SpeciesGroupSemNode: 'sp',
                          GenusGroupSemNode: 'gen',
                          VerbatimSemNode: 'verbatim',
                          CombinationSemNode: 'combin',
                          }


def _register_new_node(graph, id_minting_context, obj):
    """Returns a  canonical_id for a new obj.

    The goal is to generate unique IDs that are somewhat human readable to make it
        easier to browse the graph.
     `id_minting_context`: has
        * "parent_id" or "context_id"
        * "class_tag" or func will use the class of `obj` to tag classes
    """
    assert isinstance(id_minting_context, dict)
    pref_str, context_id = id_minting_context.get('parent_id'), ''
    if pref_str is None:
        pref_str = graph.res.base_resource.id
        context_id = id_minting_context['context_id']
    ct = id_minting_context.get('class_tag')
    if ct is None:
        ct = _NODE_CLASS_NAME_TO_CT[obj.__class__]
    can_id = _canonicalize(pref_str, ct, context_id)
    rci, n = can_id, 1
    while True:
        wtid = graph._by_id.get(can_id)
        if wtid is None:
            graph._by_id[can_id] = obj
            return can_id
        if wtid == obj:
            return can_id
        n += 1
        can_id = '{}:v{}'.format(rci, n)


def _canonicalize(res_id, pred_id, entity_id):
    ne = [str(i) for i in (res_id, pred_id, entity_id) if i]
    return ':'.join(ne)
