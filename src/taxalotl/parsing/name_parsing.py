#!/usr/bin/env python
from __future__ import print_function
import json
import logging

from .taxonomic_ranks import (SPECIES_SORTING_NUMBER,
                                      GENUS_SORTING_NUMBER)
from .parsing.gnparser import parse_name_to_dict

_LOG = logging.getLogger(__name__)


def _add_apparent_rank_and_clade_name(name_dict):
    if 'infra_epithet' in name_dict:
        name_dict['apparent_rank'] = 'infraspecies'
    elif 'sp_epithet' in name_dict:
        name_dict['apparent_rank'] = 'species'
    else:
        name_dict['apparent_rank'] = 'clade'
        if 'subgenus' in name_dict:
            name_dict['clade_name'] = name_dict['subgenus']
        elif 'genus' in name_dict:
            name_dict['clade_name'] = name_dict['genus']
        else:
            name_dict['clade_name'] = name_dict['higher_group_name']


def _convert_gnparser_detail(deet, tnd, gnp_dict):
    annot = deet.get('annotationIdentification', '')
    had_cf = False
    if annot:
        assert annot in ['sp.', 'cf.'], 'annot not sp. or cf.'
        tnd['undescribed'] = True
        if annot.startswith('cf'):
            had_cf = True
    if 'uninomial' in deet:
        tnd['higher_group_name'] = deet['uninomial']['value']
    else:
        if 'infraspecificEpithets' in deet:
            ielist = deet['infraspecificEpithets']
            assert len(ielist) == 1
            tnd['infra_epithet'] = ielist[0]['value']
        if 'specificEpithet' in deet:
            tnd['sp_epithet'] = deet['specificEpithet']['value']
        elif annot:
            tnd['sp_epithet'] = annot
        else:
            assert False, 'no specificEpithet or sp. annotationIdentification in detail'
        assert 'genus' in deet, 'genus not in detail'
        tnd['genus'] = deet['genus']['value']

    ignored = deet.get('ignored')
    if not ignored:
        ignored = gnp_dict.get('unparsedTail')
        if ignored:
            ignored = {'value': ignored.strip()}

    raw_ignored = ''
    if ignored:
        raw_ignored = ignored['value']
        tnd['specimen_code'] = raw_ignored.strip()
    if annot:
        if had_cf:
            normalized = gnp_dict.get('normalized')
        elif raw_ignored:
            normalized = gnp_dict.get('verbatim')[:-len(raw_ignored)]
        else:
            normalized = ''
        if normalized:
            tnd['normalized'] = normalized


def parse_name(name):
    """Returns name dict with the following keys:
        'undescribed' -> Absent=False. True means had "sp." or "cf "
        'combination' -> None if not below genus level or canonical name
        some of: 'higher_group_name', 'genus', 'subgenus', 'sp_epithet', 'infra_epithet',
        'specimen_code' unparsed, stripped cruft

    'apparent_rank' = 'clade' | 'species' | 'infraspecies'
        If "clade" then 'clade_name' will hold first of subgenus, genus or higher_group_name
    """
    gnp_dict = parse_name_to_dict(name)
    taxalotl_name_dict = {}
    tnd = taxalotl_name_dict
    if not gnp_dict['parsed']:
        return tnd
    emsg = 'failing on gnp_dict = {}\n from {}\n{}'
    try:
        deets = gnp_dict['details']
        if len(deets) != 1:
            assert gnp_dict['hybrid'], 'multiple details but not hybrid'
            tnd['hybrid'] = True
            return tnd
        else:
            for deet in deets:
                _convert_gnparser_detail(deet, tnd, gnp_dict)
    except Exception as x:
        raise RuntimeError(emsg.format(json.dumps(gnp_dict, indent=2), name, str(x)))
    _add_apparent_rank_and_clade_name(tnd)
    if tnd['apparent_rank'] != 'clade':
        tnd['combination'] = gnp_dict['canonicalName']['full']
        assert ' sp ' not in tnd['combination'], " sp in name"
        assert ' cf' not in tnd['combination'], " cf in name"
    return tnd


def parse_sp_name(name, rank):
    x = parse_name(name)
    return x


def parse_genus_group_name(name, rank):
    x = parse_name(name)
    if rank:
        if rank == 'genus':
            x['genus'] = x['higher_group_name']
            del x['higher_group_name']
        else:
            assert False
    return x


def parse_higher_name(name, rank):
    x = parse_name(name)
    return x


def parse_name_using_rank_hints(name, rank=None, rank_sorting_number=None):
    name_dict = None
    if rank_sorting_number is not None:
        if rank_sorting_number <= SPECIES_SORTING_NUMBER:
            name_dict = parse_sp_name(name, rank)
        elif rank_sorting_number <= GENUS_SORTING_NUMBER:
            name_dict = parse_genus_group_name(name, rank)
    if name_dict is None:
        name_dict = parse_higher_name(name, rank)

    return name_dict
