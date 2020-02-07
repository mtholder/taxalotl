#!/usr/bin/env python
from __future__ import print_function
import json
from peyotl import (get_logger, )
from .taxon_concept_node import TaxonConceptSemNode
from .verbatim_name import VerbatimSemNode
from .names_for_ranks import (GenusGroupSemNode,
                              HigherGroupSemNode,
                              SpecimenCodeSemNode,
                              SpeciesGroupSemNode,
                              TypeSpecimen
                              )
from .name import CombinationSemNode
from .graph_node import AuthoritySemNode

_LOG = get_logger(__name__)


class SemGraph(object):
    att_list = ['_authorities',
                '_combinations',
                '_genus_group_names',
                '_higher_group_names',
                '_references',
                '_species_group_epithets',
                '_specimen_codes',
                '_specimens',
                '_taxon_concepts',
                '_type_specimens',
                '_verbatim_name',
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
        self._type_specimens = None

    def impute_type_specimens(self):
        for tc in self.taxon_concept_list:
            if not tc.is_specimen_based:
                continue
            if tc.hybrid or tc.undescribed or not tc.rank:
                continue
            if tc.rank == 'species':
                epithet = tc.most_terminal_name
                try:
                    if not epithet.type_materials:
                        self._add_type_specimen(None, epithet, tc.is_synonym_of)
                except:
                    _LOG.exception("problem adding type materials")
        for tc in self.taxon_concept_list:
            if tc.hybrid or tc.undescribed or not tc.rank:
                continue
            if tc.is_specimen_based and tc.rank != 'species':
                infra_epithet = tc.most_terminal_name
                if infra_epithet is tc.has_name.sp_epithet:
                    continue
                if not infra_epithet.type_materials:
                    self._add_type_specimen(None, infra_epithet, tc.is_synonym_of)

    def denormalize_homonyms(self):
        # for species ranK:
        #   multiple authority entries
        #   same valid epithet in multiple valid genera
        dup_auth_to_mint = {}
        for tc in self.taxon_concept_list:
            if tc.rank and tc.rank == 'species':
                epithet = tc.most_terminal_name
                if epithet is None:
                    if not (tc.hybrid or tc.undescribed):
                        _LOG.warn('NO Epithet for  = {}'.format(tc.__dict__))
                    continue
                if isinstance(epithet._authority, list):
                    dup_auth_to_mint.setdefault(epithet, []).append(tc)

        for name, tc_list in dup_auth_to_mint.items():
            verb_name = tc_list[0].has_name
            other_name_tc_pairs = []
            same_name_tc = []
            for tc in tc_list[1:]:
                if tc.has_name is verb_name:
                    same_name_tc.append(tc)
                else:
                    other_name_tc_pairs.append(tc)
            for other in other_name_tc_pairs:
                self._split_tc_with_shared_sp_epithet(tc_list[0], other)
            for other in same_name_tc:
                self._split_tc_with_shared_name(tc_list[0], other)

        if self.res.id.startswith('cof'):
            import sys
            # sys.exit(1)

    def _split_tc_with_shared_sp_epithet(self, fixed, other):
        assert fixed is not other
        fix_name, oth_name = fixed.has_name, other.has_name
        _LOG.debug('splitting "{}" from ["{}"]'.format(fix_name.name, oth_name.name))
        fix_genus, oth_genus = fix_name.genus_name, oth_name.genus_name
        fix_sp_epi, oth_sp_epi = fix_name.sp_epithet, oth_name.sp_epithet
        assert fix_sp_epi is oth_sp_epi
        if fix_sp_epi in oth_genus.contained:
            oth_genus.contained.remove(fix_sp_epi)
        new_epi = self._add_sp_epithet(other, fix_sp_epi._name, oth_genus, avoid_dup=False)
        oth_genus.contained.append(new_epi)
        oth_name.sp_epithet = new_epi
        vtc = other
        if vtc._is_synonym_of:
            vtc = vtc._is_synonym_of
        for a in fix_sp_epi._authority:
            if other in a.taxon_concept_set or vtc in a.taxon_concept_set:
                new_epi.claim_authority(a)
                break
        assert new_epi._authority
        fix_sp_epi._authority.remove(new_epi._authority)
        if len(fix_sp_epi._authority) == 1:
            fix_sp_epi._authority = fix_sp_epi._authority[0]

    def _split_tc_with_shared_name(self, fixed, other):
        fix_vname = fixed.has_name
        assert fix_vname is other.has_name
        new_vname = self._add_verbatim_name(other, fix_vname.name, avoid_dup=False)
        assert not fix_vname.specimen_codes
        for attr in VerbatimSemNode.extra_pred:
            v = getattr(fix_vname, attr, None)
            if v:
                setattr(new_vname, attr, v)
        other.has_name = new_vname
        self._split_tc_with_shared_sp_epithet(fixed, other)

    def postorder_taxon_concepts(self):
        yielded = set()
        tcl = self.taxon_concept_list
        todo = []
        for tc in tcl:
            if not tc.child_set:
                yielded.add(tc)
                yield tc
            else:
                todo.append(tc)
        prev_todo_len = 1 + len(todo)
        while todo:
            assert prev_todo_len > len(todo)
            prev_todo_len = len(todo)
            ntd = []
            for tc in todo:
                if tc.child_set.issubset(yielded):
                    yielded.add(tc)
                    yield tc
                else:
                    ntd.append(tc)
            todo = ntd


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
    def canonical_name_str_to_taxon_concept_map(self):
        return {i.canonical_name.name: i for i in self._taxon_concepts}

    @property
    def valid_name_to_taxon_concept_map(self):
        return {i.valid_name.name: i for i in self.valid_taxa_dict.values()}

    def get_by_id(self, can_id, default=None):
        return self._by_id.get(can_id, default)

    def _add_name(self, container, node_type, parent_sem_node, name, extra_container=None, avoid_dup=True):
        search_cont = container if extra_container is None else extra_container.contained
        x = None if (not avoid_dup) else _find_by_name(search_cont, name)
        if x is None:
            d = {'parent_id': parent_sem_node.canonical_id}
            if extra_container is not None:
                d['class_tag'] = 'epi'  # need to figure out if this is the best choice for any extra container obj
            x = node_type(self, d, name)
            search_cont.append(x)
            if search_cont is not container:
                container.append(x)
        return x

    def add_authority(self, tax_con_sem_node, name_sem, authors, year):
        auth_list = self.authorities
        x = None
        for a in auth_list:
            if a.authors == authors and a.year == year:
                x = a
                break
        if x is None:
            d = {'parent_id': name_sem.canonical_id}
            x = AuthoritySemNode(self, d, authors, year, tax_con_sem_node)
        else:
            x.taxon_concept_set.add(tax_con_sem_node)
        auth_list.append(x)
        name_sem.claim_authority(x)
        return x

    def _add_normalized(self, par_sem_node, name_str):
        return self._add_name(self.combinations, CombinationSemNode, par_sem_node, name_str)

    def _add_combination(self, par_sem_node, name_str):
        return self._add_name(self.combinations, CombinationSemNode, par_sem_node, name_str)

    def _add_verbatim_name(self, tax_con_sem_node, name_str, avoid_dup=True):
        return self._add_name(self.verbatim_name, VerbatimSemNode, tax_con_sem_node, name_str, avoid_dup=avoid_dup)

    def _add_genus(self, par_sem_node, name_str):
        return self._add_name(self.genus_group_names, GenusGroupSemNode, par_sem_node, name_str)

    _add_subgenus = _add_genus

    def _add_higher_group_name(self, par_sem_node, name_str):
        return self._add_name(self.higher_group_names, HigherGroupSemNode, par_sem_node, name_str)

    def _add_sp_epithet(self, par_sem_node, name_str, prev_word_sn, avoid_dup=True):
        return self._add_name(self.species_group_epithets, SpeciesGroupSemNode,
                              par_sem_node, name_str, prev_word_sn, avoid_dup=avoid_dup)

    _add_infra_epithet = _add_sp_epithet

    def _add_specimen_code(self, par_sem_node, name_str):
        return self._add_name(self.specimen_codes, SpecimenCodeSemNode, par_sem_node, name_str)

    def add_taxon_concept(self, foreign_id):
        x = TaxonConceptSemNode(self, foreign_id)
        self.taxon_concepts.append(x)
        return x

    def remove_taxon_concept(self, tc):
        if tc in self._taxon_concepts:
            self._taxon_concepts.remove(tc)
        if tc.canonical_id in self._by_id:
            del self._by_id[tc.canonical_id]

    def _add_type_specimen(self, spec_code, epithet_syn_name, valid_taxon):
        d = {'parent_id': epithet_syn_name.canonical_id}
        x = TypeSpecimen(self, d, spec_code, epithet_syn_name, valid_taxon)
        self.type_specimens.append(x)
        epithet_syn_name.claim_type_material(x)
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
