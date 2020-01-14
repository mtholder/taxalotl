#!/usr/bin/env python
from __future__ import print_function

import re
from peyotl import get_logger

_LOG = get_logger(__name__)

_SP_EPITHET = re.compile(r'^([A-Z][a-z]+)\s+sp[.]?\s*$')
_CF_SPECIES = re.compile(r'^([A-Z][a-z]+)\s+cf[.]?\s*([a-z]+)$')
_BINOMINAL = re.compile(r'^([A-Z][a-z]+)\s+([a-z]+)\s*$')
_SP_EPITHET_CRUFT = re.compile(r'^([A-Z][a-z]+)\s+sp[.]?\s+(\S.*)$')
_CF_SPECIES_CRUFT = re.compile(r'^([A-Z][a-z]+)\s+cf[.]?\s*([a-z]+)\s+(\S.*)$')
_BINOMINAL_CRUFT = re.compile(r'^([A-Z][a-z]+)\s+([a-z]+)\s+(\S.*)$')
_UNINOMINAL = re.compile(r'^([A-Z][a-z]+)$')
_TWO_CAP = re.compile(r'^([A-Z][a-z]+)\s+([A-Z][a-z]+)$')
_PARENS_SUB_GENUS = re.compile(r'^([A-Z][a-z]+)\s+\(([A-Z][a-z]+)\)$')


def _check_genus_sp_dot(name, d):
    m = _SP_EPITHET.match(name)
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
    m = _SP_EPITHET_CRUFT.match(name)
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
    d['undescribed'] = m.group(3).strip()
    return True


_TRINOMINAL = re.compile(r'^([A-Z][a-z]+)\s+([a-z]+)\s+([a-z]+)\s*$')


def _check_trinomial(name, d):
    m = _TRINOMINAL.match(name)
    if not m:
        return False
    d['genus'] = m.group(1)
    d['sp_epithet'] = m.group(2)
    d['infra_epithet'] = m.group(3)
    return True


sp_check_order = (_check_genus_sp_dot, _check_cf_species, _check_genus_sp_epithet,
                  _check_genus_sp_dot_cruft, _check_cf_species_cruft, _check_genus_sp_epithet_cruft)

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
