#!/usr/bin/env python
from __future__ import print_function

import json
import os

from peyotl import (get_logger, write_as_json)

from taxalotl.util import OutFile, OutDir
from taxalotl.taxonomic_ranks import (ABOVE_GENUS_SORTING_NUMBER,
                                      SPECIES_SORTING_NUMBER,
                                      GENUS_SORTING_NUMBER,
                                      MINIMUM_HIGHER_TAXON_NUMBER)

from taxalotl.parsing.name_parsing import parse_name_using_rank_hints, parse_name_to_dict
from taxalotl.tax_partition import IGNORE_COMMON_NAME_SYN_TYPES
from taxalotl.sem_graph.graph import SemGraph

_LOG = get_logger(__name__)


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
            sem_node.claim_uninomial_synonym(syn.name, syn_type=syn.syn_type, **pnd)
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
    for syn in valid_taxon.synonym_list:
        shn = syn.has_name
        if shn and shn.name == fn:
            return shn
    raise RuntimeError('"{}" is not a basionym'.format(fn))


def add_authority_to_name(tax_con_sem_node, name_sem, authors, year):
    mtn = name_sem.most_terminal_name
    # _LOG.debug('auth for {} is {}'.format(mtn, authors))
    tax_con_sem_node.graph.add_authority(tax_con_sem_node, mtn, authors, year)


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
    add_authority_to_name(taxon_sem_node, name_sem, authors, year)


def semanticize_node_auth_synonym(res, sem_graph, node, sem_node, syn):
    shn = sem_node.has_name
    _parse_auth_syn(res, sem_graph, sem_node, syn)


def semanticize_node_name(res, sem_graph, tc, node):
    try:
        rsn = node.rank_sorting_number()
    except KeyError:
        rsn = None
    name_dict = parse_name_using_rank_hints(node.name, node.rank, rsn)
    semanticize_names(tc, node.name, name_dict, node)


def semanticize_names(tax_con_sem_node, name, name_dict, node=None):
    """Translates content in the parsed `name_dict` in to graph nodes for this taxon_concept_sem_node
    If `node` is not None, it should be a `Taxon` instance (used for flags)
    """
    tcsn = tax_con_sem_node
    undescribed = name_dict.get('undescribed')
    if undescribed:
        tcsn.claim_undescribed()
    if name_dict.get('hybrid'):
        tcsn.claim_hybrid()
    vnsn = tcsn.add_verbatim_name(name)
    # _LOG.debug('semanticizing {} for {}'.format(name, tcsn.canonical_id))
    valid_tcsn = tcsn if tcsn._is_synonym_of is None else tcsn._is_synonym_of
    valid_nph = valid_tcsn.has_name
    if undescribed:
        normalized = name_dict.get('normalized')
        if normalized:
            vnsn.add_normalized(normalized)
    else:
        combination = name_dict.get('combination')
        if combination:
            vnsn.add_combination(combination)
    genus = name_dict.get('genus')
    if genus:
        vnsn.add_genus(genus)
    valid_genus = valid_nph.genus_name
    subgenus = name_dict.get('subgenus')
    if subgenus:
        vnsn.add_subgenus(subgenus)
    sp_epithet = name_dict.get('sp_epithet')
    valid_sp = None
    if sp_epithet:
        assert valid_genus is not None
        valid_sp = vnsn.add_sp_epithet(sp_epithet, valid_genus)
    infra_epithet = name_dict.get('infra_epithet')
    if infra_epithet:
        assert valid_sp is not None
        vnsn.add_infra_epithet(infra_epithet, valid_sp)
    higher_group_name = name_dict.get('higher_group_name')
    if higher_group_name:
        vnsn.add_higher_group_name(higher_group_name)
    specimen_code = name_dict.get('specimen_code')
    if specimen_code:
        vnsn.add_specimen_code(specimen_code)
    if node and node.flags:
        for flag in node.flags:
            if flag == 'infraspecific':
                if infra_epithet is None:
                    m = 'taxon "{}" flagged as infraspecific, but not parsed as such'
                    raise ValueError(m.format(name))
            else:
                tcsn.claim_flag(flag)
    return vnsn


def semanticize_and_serialize_tax_part(taxolotl_config, res, fragment, out_dir, tax_part, tax_forest):
    sem_graph = semanticize_tax_part(taxolotl_config, res, fragment, tax_part, tax_forest)
    serialize_sem_graph(taxolotl_config, sem_graph, out_dir)
    return sem_graph


# noinspection PyUnusedLocal
def semanticize_tax_part(taxolotl_config, res, fragment, tax_part, tax_forest):
    sem_graph = SemGraph(taxolotl_config, res)
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
    sem_graph.denormalize_homonyms()
    sem_graph.impute_type_specimens()
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
