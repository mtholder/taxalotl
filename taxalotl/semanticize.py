#!/usr/bin/env python
from __future__ import print_function

import json
import os

from peyotl import (get_logger, write_as_json)

from taxalotl.util import OutFile, OutDir
from .taxonomic_ranks import (ABOVE_GENUS_SORTING_NUMBER,
                              SPECIES_SORTING_NUMBER,
                              GENUS_SORTING_NUMBER,
                              MINIMUM_HIGHER_TAXON_NUMBER)

from .name_parsing import parse_name_using_rank_hints, parse_name_to_dict
from .tax_partition import IGNORE_COMMON_NAME_SYN_TYPES


_LOG = get_logger(__name__)


def serialize_triple_object(o):
    return o.canonical_id if isinstance(o, SemGraphNode) else o


class SemGraphNode(object):
    def __init__(self, sem_graph, canoncial_id):
        self.canonical_id = sem_graph.register_obj(canoncial_id, self)
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
    if pred_id:
        return '{}:{}:{}'.format(res_id, pred_id, entity_id)
    return '{}:{}'.format(res_id, pred_id, entity_id)

KNOWN_FLAGS = frozenset(['hidden', 'sibling_higher'])

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
        self._is_synonym_of = None
        self.problematic_synonyms = None
        self.synonyms = None
        self.syn_type = None
        self.former_ranks = None
        self.hybrid = None
        self.incertae_sedis = None
        self.other_flags = None

    def explain(self, out):
        cn = self.canonical_name.name if self.canonical_name else ''
        out.write('"{}" Taxon({}) '.format(cn, self.canonical_id))
        if self._is_synonym_of:
            out.write('synonym of "{}"'.format(self.valid_name.name))
        else:
            out.write(' rank={}'.format(self.rank if self.rank else '?'))


    @property
    def valid_name(self):
        if self._is_synonym_of:
            return self._is_synonym_of.valid_name
        return self.canonical_name

    @property
    def canonical_name(self):
        if not self.has_name:
            return None
        return self.has_name.canonical_name

    @property
    def is_synonym_of(self):
        return bool(self._is_synonym_of)

    @property
    def is_specimen_based(self):
        return self.rank in ['species', 'subspecies']

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

    def claim_hybrid(self):
        self.hybrid = True

    def claim_incertae_sedis(self):
        self.incertae_sedis = True

    def claim_flag(self, flag):
        if flag == 'incertae_sedis':
            self.claim_incertae_sedis
        elif flag == 'hybrid':
            self.claim_hybrid()
        else:
            if flag not in KNOWN_FLAGS:
                raise ValueError('Unknown flag "{}"'.format(flag))
            self._add_other_flag(flag)

    def _add_other_flag(self, flag):
        if self.other_flags:
            if flag in self.other_flags:
                return
            self.other_flags.append(flag)
            self.other_flags.sort()
        else:
            self.other_flags = [flag]

    @property
    def predicates(self):
        return ['hybrid', 'is_child_of', 'rank', 'has_name', 'id', 'undescribed',
                'is_synonym',
                'incertae_sedis', 'other_flags'
                'syn_type', 'former_ranks',
                'problematic_synonyms', 'synonyms']

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

    def _get_next_syn_id(self):
        ns = len(self.synonyms) if self.synonyms else 0
        return '{}:syn{}'.format(self.concept_id, ns)

    @property
    def is_synonym(self):
        return self._is_synonym_of is not None

    def claim_is_synonym_of(self, valid):
        self._is_synonym_of = valid

    def claim_former_rank(self, rank):
        if self.former_ranks is None:
            self.former_ranks = []
        if rank not in self.former_ranks:
            self.former_ranks.append(rank)

    def claim_uninomial_synonym(self, res, name, syn_type, **kwargs):
        if 'genus' in kwargs:
            return self._claim_syn_impl(res, name, syn_type, 'genus', **kwargs)
        if 'higher_group_name' not in kwargs:
            raise ValueError("Expecting 'higher_group_name' or 'genus' kwarg "
                             "in claim_uninomial_synonym")
        return self._claim_syn_impl(res, name, syn_type, 'clade', **kwargs)

    def claim_formerly_full_species(self, res, name, syn_type, **kwargs):
        self.claim_former_rank('species')
        return self.claim_binom_synonym(res, name, syn_type, **kwargs)

    def claim_formerly_subspecies(self, res, name, syn_type, **kwargs):
        self.claim_former_rank('subspecies')
        return self.claim_trinomial_synonym(res, name, syn_type, **kwargs)

    def claim_binom_synonym(self, res, name, syn_type, **kwargs):
        for expected in ('genus', 'sp_epithet'):
            if expected not in kwargs:
                raise ValueError("Expecting '{}' kwarg in claim_binom_synonym".format(expected))
        return self._claim_syn_impl(res, name, syn_type, 'species', **kwargs)

    def claim_trinomial_synonym(self, res, name, syn_type, **kwargs):
        for expected in ('genus', 'sp_epithet', 'infra_epithet'):
            if expected not in kwargs:
                raise ValueError("Expecting '{}' kwarg in claim_trinomial_synonym".format(expected))
        return self._claim_syn_impl(res, name, syn_type, 'infraspecies', **kwargs)

    def _add_to_syn_list(self, syntc):
        if self.synonyms is None:
            self.synonyms = []
        self.synonyms.append(syntc)

    def _claim_syn_impl(self, res, name, syn_type, rank, **kwargs):
        tc = self.graph.add_taxon_concept(res, self._get_next_syn_id())
        self._add_to_syn_list(tc)
        tc.claim_rank(rank)
        tc.claim_is_synonym_of(self)
        if syn_type:
            tc.syn_type = syn_type
        if kwargs.get('undescribed', False):
            tc.claim_undescribed()
        semanticize_names(res, self.graph, tc, name, kwargs, None)
        return tc

class AuthoritySemNode(SemGraphNode):
    auth_sem_nd_pred = ('authors', 'year')

    def __init__(self, sem_graph, res_id, tag, authors, year):
        concept_id = 'auth'
        ci = canonicalize(tag, '', concept_id)
        super(AuthoritySemNode, self).__init__(sem_graph, ci)
        self._authors = authors
        self._year = year

    @property
    def year(self):
        return self._year

    @property
    def authors(self):
        return self._authors

    @property
    def predicates(self):
        return AuthoritySemNode.auth_sem_nd_pred

class NameSemNode(SemGraphNode):
    name_sem_nd_pred = ('name', 'authority')

    def __init__(self, sem_graph, res_id, tag, concept_id, name):
        ci = canonicalize(res_id, tag, concept_id)
        super(NameSemNode, self).__init__(sem_graph, ci)
        self._name = name
        self._authority = None

    @property
    def name(self):
        return self._name

    @property
    def predicates(self):
        return NameSemNode.name_sem_nd_pred

    def claim_authority(self, auth):
        self._authority = None

    @property
    def authority(self):
        return None if self._authority is None else self._authority.canonical_id


class CombinationSemNode(NameSemNode):
    def __init__(self, sem_graph, res_id, concept_id, name):
        super(CombinationSemNode, self).__init__(sem_graph, res_id, 'combin', concept_id, name)


class VerbatimSemNode(NameSemNode):
    extra_pred = ('combination', 'higher_group_name', 'genus_name', 'authority',
                  'subgenus_names', 'sp_epithet', 'infra_epithets', 'specimen_codes')
    name_sem_nd_pred = tuple(list(NameSemNode.name_sem_nd_pred) + list(extra_pred))

    def __init__(self, sem_graph, res_id, concept_id, name):
        super(VerbatimSemNode, self).__init__(sem_graph, res_id, 'verbatim', concept_id, name)
        self.combination = None
        self.higher_group_name = None
        self.genus_name = None
        self.subgenus_names = None
        self.sp_epithet = None
        self.infra_epithets = None
        self.specimen_codes = None

    @property
    def canonical_name(self):
        if self.combination is not None:
            return self.combination
        if self.subgenus_names is not None:
            return self.subgenus_names
        if self.genus_name is not None:
            return self.genus_name
        if self.higher_group_name is not None:
            return self.higher_group_name
        assert self.name
        return self


    @property
    def predicates(self):
        return VerbatimSemNode.name_sem_nd_pred

    def claim_higher_group_name(self, n):
        assert self.higher_group_name is None
        self.higher_group_name = n

    def claim_combination(self, n):
        assert self.combination is None or self.combination is n
        self.combination = n

    def claim_genus(self, n):
        assert self.genus_name is None or self.genus_name == n
        self.genus_name = n

    def claim_subgenus(self, n):
        if self.subgenus_names is None:
            self.subgenus_names = [n]
        else:
            self.subgenus_names.append(n)

    def claim_sp_epithet(self, n):
        assert self.sp_epithet is None or self.sp_epithet is n
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
    def most_terminal_name(self):
        ie = self.most_terminal_infra_epithet
        if ie is not None:
            return ie
        for n in [self.sp_epithet, self.subgenus_names, self.genus_name, self.higher_group_name]:
            if n:
                return n
        return None

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
        self.contained = []


class SpeciesGroupSemNode(NameSemNode):
    sp_grp_name_sem_nd_pred = tuple(list(NameSemNode.name_sem_nd_pred) + ['type_materials', 'authority'])

    def __init__(self, sem_graph, res_id, concept_id, name):
        super(SpeciesGroupSemNode, self).__init__(sem_graph, res_id, 'sp', concept_id, name)
        self.type_materials = None
        self._authority = None

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

    def __init__(self):
        self._by_id = {}
        for att in SemGraph.att_list:
            setattr(self, att, None)

    def _all_specimen_based_tax_con_dict(self):
        raw = self._taxon_concepts if self._taxon_concepts else []
        r = {}
        for tcobj in raw:
            if tcobj.is_specimen_based:
                r[tcobj.canonical_id] = tcobj
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

    def register_obj(self, can_id, obj):
        rci, n = can_id, 1
        while True:
            wtid = self._by_id.get(can_id)
            if wtid is None:
                self._by_id[can_id] = obj
                return can_id
            if wtid == obj:
                return can_id
            n += 1
            can_id = '{}:v{}'.format(rci, n)

    def get_by_id(self, can_id, default=None):
        return self._by_id.get(can_id, default)

    def _add_name(self, container, node_type, res_id, concept_id, name, extra_container=None):
        if extra_container is None:
            search_cont, id_minting_str = container, concept_id
        else:
            search_cont = extra_container.contained
            pref = '{}:'.format(res_id)
            assert extra_container.canonical_id.startswith(pref)
            sans_res = extra_container.canonical_id[len(pref):]
            id_minting_str = '{}.epi.{}'.format(sans_res, concept_id)
        x = _find_by_name(search_cont, name)
        if x is None:
            x = node_type(self, res_id, id_minting_str, name)
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
            x = AuthoritySemNode(self, res_id, name_sem.canonical_id, authors, year)
        auth_list.append(x)
        name_sem.claim_authority(x)
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

    def add_sp_epithet(self, res_id, concept_id, name, genus):
        return self._add_name(self.species_group_epithets, SpeciesGroupSemNode,
                              res_id, concept_id, name, genus)

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
    pnd = parse_name_using_rank_hints(syn.name)
    if not pnd:
        sem_node.claim_problematic_synonym_statement(syn.name, syn.syn_type, "unparseable")
        return
    arank = pnd['apparent_rank']
    if arank == 'clade':
        if rsn <= SPECIES_SORTING_NUMBER:
            e = "higher taxon name synonym for taxon rank = {}".format(node.rank)
            sem_node.claim_problematic_synonym_statement(syn.name, syn.syn_type, e)
        else:
            cn = pnd['clade_name']
            del pnd['clade_name']
            if rsn == GENUS_SORTING_NUMBER:
                pnd['genus'] = cn
            else:
                pnd['higher_group_name'] = cn
            sem_node.claim_uninomial_synonym(res, syn.name, syn_type=syn.syn_type, **pnd)
    elif arank == 'species':
        if rsn > SPECIES_SORTING_NUMBER:
            e = "species synonym for taxon rank = {}".format(node.rank)
            sem_node.claim_problematic_synonym_statement(syn.name, syn.syn_type, e)
        elif rsn == SPECIES_SORTING_NUMBER:
            sem_node.claim_binom_synonym(res, syn.name, syn_type=syn.syn_type, **pnd)
            # ['genus'], pnd['sp_epithet'], pnd.get('undescribed', False),
            #                             syn.syn_type, pnd.get('specimen_code'))
        else:
            sem_node.claim_formerly_full_species(res, syn.name, syn_type=syn.syn_type, **pnd)
    else:
        assert arank == 'infraspecies'
        if rsn > SPECIES_SORTING_NUMBER:
            e = "infraspecies synonym for taxon rank = {}".format(node.rank)
            sem_node.claim_problematic_synonym_statement(syn.name, syn.syn_type, e)
        elif rsn == SPECIES_SORTING_NUMBER:
            sem_node.claim_formerly_subspecies(res, syn.name, syn_type=syn.syn_type, **pnd)
        else:
            sem_node.claim_trinomial_synonym(res, syn.name, syn_type=syn.syn_type, **pnd)
    '''

    

    else:
        valid_combo = sem_node.valid_combination
        if len(syn.name.split()) != len(valid_combo.split()):
            _LOG.debug('WARNING: "{}" is a {} for {} ({})'.format(syn.name, syn.syn_type, node.name, node.id))
        else:
            _LOG.debug('         "{}" is a {} for {} ({})'.format(syn.name, syn.syn_type, node.name, node.id))
    st = syn.syn_type
    '''

def _assure_nonbasionym_exists_as_syn(res, sem_graph, valid_taxon, canonical_name):
    fn = canonical_name['full']
    vthn = valid_taxon.has_name
    if vthn and vthn.canonical_name and vthn.canonical_name.name == fn:
        return vthn
    if valid_taxon.synonyms:
        for syn in valid_taxon.synonyms:
            shn = syn.has_name
            if shn and shn.name == fn:
                return shn
    raise RuntimeError('"{}" is not a basionym'.format(fn))

def add_authority_to_name(res, sem_graph, name_sem, authors, year):
    mtn = name_sem.most_terminal_name
    # _LOG.debug('auth for {} is {}'.format(mtn, authors))
    sem_graph.add_authority(res.id, mtn, authors, year)


def _assure_basionym_exists_as_syn(res, sem_graph, valid_taxon, canonical_name):
    fn = canonical_name['full']
    vthn = valid_taxon.has_name
    if vthn and vthn.canonical_name and vthn.canonical_name.name == fn:
        return vthn
    raise RuntimeError('"{}" is a basionym'.format(fn))

def _parse_auth_syn(res, sem_graph, taxon_sem_node, syn):
    gnp = parse_name_to_dict(syn.name)
    assert gnp['parsed']
    cn = gnp['canonicalName']
    details = gnp['details']
    try:
        if isinstance(details, list):
            assert len(details) == 1
            details = details[0]
        if 'infraspecificEpithets' in details:
            target_name = details['infraspecificEpithets'][-1]
        elif 'specificEpithet' in details:
            target_name = details['specificEpithet']
        else:
            assert 'uninomial' in details
            target_name = details['uninomial']
    except:
        import sys
        sys.exit('choking on\n{}\n'.format(json.dumps(details, indent=2)))
    if isinstance(target_name, list):
        assert len(target_name) == 1
        target_name = target_name[0]
    authorship = target_name['authorship']
    parens_preserved = authorship['value'].strip()
    is_basionym = parens_preserved.startswith('(')
    if is_basionym:
        name_sem = _assure_basionym_exists_as_syn(res, sem_graph, taxon_sem_node, cn)
    else:
        name_sem = _assure_nonbasionym_exists_as_syn(res, sem_graph, taxon_sem_node, cn)
    ba = authorship.get('basionymAuthorship', {})
    assert ba
    authors, year = ba.get('authors', []), ba.get('year', {}).get('value')
    add_authority_to_name(res, sem_graph, name_sem, authors, year)


def semanticize_node_auth_synonym(res, sem_graph, node, sem_node, syn):
    shn = sem_node.has_name
    _parse_auth_syn(res, sem_graph, sem_node, syn)

def semanticize_node_name(res, sem_graph, tc, node):
    try:
        rsn = node.rank_sorting_number()
    except KeyError:
        rsn = None
    name_dict = parse_name_using_rank_hints(node.name, node.rank, rsn)
    semanticize_names(res, sem_graph, tc, node.name, name_dict, node)


def semanticize_names(res, sem_graph, taxon_concept_sem_node, name, name_dict, node=None):
    """
    """
    tcsn = taxon_concept_sem_node
    if name_dict.get('undescribed'):
        tcsn.claim_undescribed()
    if name_dict.get('hybrid'):
        tcsn.claim_hybrid()


    rn = sem_graph.add_verbatim_name(res.base_resource.id, tcsn.concept_id, name)
    _LOG.debug('semanticizing {} for {}'.format(name, tcsn.concept_id))
    name_part_holder = rn
    tcsn.claim_name(rn)
    valid_tcsn = tcsn if tcsn._is_synonym_of is None else tcsn._is_synonym_of
    valid_nph = valid_tcsn.has_name
    combination = name_dict.get('combination')
    bresid = res.base_resource.id
    if combination:
        cn = sem_graph.add_combination(bresid, tcsn.concept_id, combination)
        rn.claim_combination(cn)
    genus = name_dict.get('genus')
    if genus:
        cn = sem_graph.add_genus(bresid, tcsn.concept_id, genus)
        name_part_holder.claim_genus(cn)
    valid_genus = valid_nph.genus_name
    subgenus = name_dict.get('subgenus')
    if subgenus:
        cn = sem_graph.add_subgenus(bresid, tcsn.concept_id, subgenus)
        name_part_holder.claim_subgenus(cn)
    sp_epithet = name_dict.get('sp_epithet')
    if sp_epithet:
        assert valid_genus is not None
        cn = sem_graph.add_sp_epithet(bresid, tcsn.concept_id, sp_epithet, valid_genus)
        name_part_holder.claim_sp_epithet(cn)
    infra_epithet = name_dict.get('infra_epithet')
    if infra_epithet:
        assert valid_genus is not None
        cn = sem_graph.add_infra_epithet(bresid, tcsn.concept_id, infra_epithet, valid_genus)
        name_part_holder.claim_infra_epithet(cn)
    higher_group_name = name_dict.get('higher_group_name')
    if higher_group_name:
        cn = sem_graph.add_higher_group_name(bresid, tcsn.concept_id, higher_group_name)
        name_part_holder.claim_higher_group_name(cn)
    specimen_code = name_dict.get('specimen_code')
    if specimen_code:
        cn = sem_graph.add_specimen_code(bresid, tcsn.concept_id, specimen_code)
        name_part_holder.claim_specimen_code(cn)
    if node and node.flags:
        for flag in node.flags:
            if flag == 'infraspecific':
                if infra_epithet is None:
                    m = 'taxon "{}" flagged as infraspecific, but not parsed as such'
                    raise ValueError(m.format(name))
            else:
                tcsn.claim_flag(flag)
    return rn


def semanticize_and_serialize_tax_part(taxolotl_config, res, fragment, out_dir, tax_part, tax_forest):
    sem_graph = semanticize_tax_part(taxolotl_config, res, fragment, tax_part, tax_forest)
    serialize_sem_graph(taxolotl_config, sem_graph, out_dir)
    return sem_graph


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
