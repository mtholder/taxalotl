#!/usr/bin/env python
from __future__ import print_function

import re
from peyotl import get_logger

_LOG = get_logger(__name__)

_BINOMINAL = re.compile(r'^([A-Z][a-z]+)\s+([a-z]+)\s*$')
_BINOMINAL_CRUFT = re.compile(r'^([A-Z][a-z]+)\s+([a-z]+)\s+(\S.*)$')
_CF_SPECIES = re.compile(r'^([A-Z][a-z]+)\s+cf[.]?\s*([a-z]+)$')
_CF_SPECIES_CRUFT = re.compile(r'^([A-Z][a-z]+)\s+cf[.]?\s*([a-z]+)\s+(\S.*)$')
_PARENS_SUB_GENUS = re.compile(r'^([A-Z][a-z]+)\s+\(([A-Z][a-z]+)\)$')
_SPDOT_EPITHET = re.compile(r'^([A-Z][a-z]+)\s+sp[.]?\s*$')
_SPDOT_EPITHET_CRUFT = re.compile(r'^([A-Z][a-z]+)\s+sp[.]?\s+(\S.*)$')
_TRINOMINAL = re.compile(r'^([A-Z][a-z]+)\s+([a-z]+)\s+([a-z]+)\s*$')
_TRINOMINAL_CRUFT = re.compile(r'^([A-Z][a-z]+)\s+([a-z]+)\s+([a-z]+)\s+(\S.*)$')
_TWO_CAP = re.compile(r'^([A-Z][a-z]+)\s+([A-Z][a-z]+)$')
_UNINOMINAL = re.compile(r'^([A-Z][a-z]+)$')
_UNINOMINAL_CRUFT = re.compile(r'^([A-Z][a-z]+)\s+(\S.*)$')


def _check_genus_sp_dot(name, d):
    m = _SPDOT_EPITHET.match(name)
    if not m:
        return False
    d['genus'] = m.group(1)
    d['sp_epithet'] = 'sp.'
    d['undescribed'] = True
    return True


def _check_cf_species(name, d):
    m = _CF_SPECIES.match(name)
    if not m:
        return False
    d['genus'] = m.group(1)
    d['sp_epithet'] = 'cf. {}'.format(m.group(2))
    d['undescribed'] = True
    return True


def _check_genus_sp_epithet(name, d):
    m = _BINOMINAL.match(name)
    if not m:
        return False
    d['genus'] = m.group(1)
    d['sp_epithet'] = m.group(2)
    return True


def _check_genus_sp_dot_cruft(name, d):
    m = _SPDOT_EPITHET_CRUFT.match(name)
    if not m:
        return False
    d['genus'] = m.group(1)
    d['sp_epithet'] = 'sp.'
    d['undescribed'] = True
    d['specimen_code'] = m.group(2).strip()
    return True


def _check_cf_species_cruft(name, d):
    m = _CF_SPECIES_CRUFT.match(name)
    if not m:
        return False
    d['genus'] = m.group(1)
    d['sp_epithet'] = 'cf. {}'.format(m.group(2))
    d['undescribed'] = True
    d['specimen_code'] = m.group(3).strip()
    return True


def _check_genus_sp_epithet_cruft(name, d):
    m = _BINOMINAL_CRUFT.match(name)
    if not m:
        return False
    d['genus'] = m.group(1)
    d['sp_epithet'] = m.group(2)
    d['specimen_code'] = m.group(3).strip()
    return True


def _check_trinomial(name, d):
    m = _TRINOMINAL.match(name)
    if not m:
        return False
    d['genus'] = m.group(1)
    d['sp_epithet'] = m.group(2)
    d['infra_epithet'] = m.group(3)
    return True

def _check_trinomial_cruft(name, d):
    m = _TRINOMINAL.match(name)
    if not m:
        return False
    d['genus'] = m.group(1)
    d['sp_epithet'] = m.group(2)
    d['infra_epithet'] = m.group(3)
    d['specimen_code'] = m.group(4).strip()
    return True



just_sp_check_order = (_check_genus_sp_dot, _check_cf_species, _check_genus_sp_epithet)
sp_with_cruft_order = (_check_genus_sp_dot_cruft, _check_cf_species_cruft, _check_genus_sp_epithet_cruft)
sp_check_order = (_check_genus_sp_dot, _check_cf_species, _check_genus_sp_epithet,
                  _check_genus_sp_dot_cruft, _check_cf_species_cruft, _check_genus_sp_epithet_cruft)

sub_sp_check_order = (_check_trinomial,)
infra_sp_check_order = (_check_trinomial,)


def parse_sp_name(name, rank):
    d = {}
    found = False
    if rank == 'species':
        for check in sp_check_order:
            if check(name, d):
                found = True
                break
    else:
        for check in infra_sp_check_order:
            if check(name, d):
                found = True
                break
    if found:
        if rank == 'species':
            d['combination'] = '{} {}'.format(d['genus'], d['sp_epithet'])
            return d
        if d.get('infra_epithet'):
            d['combination'] = '{} {} {}'.format(d['genus'], d['sp_epithet'], d['infra_epithet'])
            return d
    raise ValueError('Could not parse species level name "{}" with rank={}'.format(name, rank))


def parse_as_just_clade_name(name, dest):
    m = _UNINOMINAL.match(name)
    if m:
        dest['apparent_rank'] = 'clade'
        dest['clade_name'] = m.group(1)
        return dest


def parse_as_just_sp_name(name, dest):
    for check in just_sp_check_order:
        if check(name, dest):
            dest['apparent_rank'] = 'species'
            return dest


def parse_as_just_subsp_name(name, dest):
    for check in sub_sp_check_order:
        if check(name, dest):
            dest['apparent_rank'] = 'infraspecies'
            return dest


def parse_as_just_infra_sp(name, dest):
    raise NotImplemented('parse_as_just_infra_sp')
    # for check in infra_sp_check_order:
    #     if check(name, dest):
    #        dest['apparent_rank'] = 'infraspecies'
    #         return dest


def parse_as_infra_sp_with_codes(name, dest):
    raise NotImplemented('parse_as_infra_sp_with_codes')
    # for check in (_TRINOMINAL_CRUFT, ):
    #     if check(name, dest):
    #        dest['apparent_rank'] = 'infraspecies'
    #         return dest


def parse_as_subsp_name_with_codes(name, dest):
    for check in (_check_trinomial_cruft,):
        if check(name, dest):
            dest['apparent_rank'] = 'infraspecies'
            return dest


def parse_as_sp_with_codes(name, dest):
    for check in sp_with_cruft_order:
        if check(name, dest):
            dest['apparent_rank'] = 'species'
            return dest


def parse_as_clade_with_codes(name, dest):
    m = _UNINOMINAL_CRUFT.match(name)
    if m:
        dest['apparent_rank'] = 'clade'
        dest['clade_name'] = m.group(1)
        dest['specimen_code'] = m.group(2).strip()
        return dest


_without_context_order = (parse_as_just_clade_name,
                          parse_as_just_sp_name,
                          parse_as_just_subsp_name,
                          # parse_as_just_infra_sp, # NotImplemented
                          # parse_as_infra_sp_with_codes,
                          parse_as_subsp_name_with_codes,
                          parse_as_sp_with_codes,
                          parse_as_clade_with_codes,)


def parse_name_string_without_context(name):
    dest = {}
    for pfn in _without_context_order:
        if pfn(name, dest):
            break
    return dest


def parse_genus_group_name(name, rank):
    d = {}
    genus, sub_genus = None, None
    if rank == 'subgenus':
        m = _TWO_CAP.match(name)
        if not m:
            m = _PARENS_SUB_GENUS.match(name)
        if m:
            genus = m.group(1)
            sub_genus = m.group(2)
    else:
        m = _UNINOMINAL.match(name)
        if m:
            genus = m.group(1)
    if genus is not None:
        d['genus'] = genus
        if rank == 'subgenus':
            if sub_genus is not None:
                d['subgenus'] = sub_genus
                d['combination'] = '{} ({})'.format(genus, sub_genus)
                return d
        else:
            return d
    raise ValueError('Could not parse genus level name "{}" with rank={}'.format(name, rank))


def parse_higher_name(name, rank):
    d = {}
    taxon_name = None
    m = _UNINOMINAL.match(name)
    if m:
        taxon_name = m.group(1)
    if taxon_name is not None:
        d['higher_group_name'] = taxon_name
        return d
    raise ValueError('Could not parse higher group name "{}" with rank={}'.format(name, rank))
