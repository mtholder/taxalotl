#!/usr/bin/env python
# -*- coding: utf-8 -*-

RANK_TO_SORTING_NUMBER = {
    "superkingdom":    2900,
    "kingdom":         2800,
    "subkingdom":      2700,
    "superphylum":     2600,
    "phylum":          2500,
    "subphylum":       2400,
    "superclass":      2300,
    "class":           2200,
    "subclass":        2100,
    "infraclass":      2000,
    "cohort":          1900,
    "superorder":      1800,
    "order":           1700,
    "suborder":        1600,
    "infraorder":      1500,
    "parvorder":       1400,
    "superfamily":     1300,
    "family":          1200,
    "subfamily":       1100,
    "tribe":           1000,
    "subtribe":         900,
    "genus":            800,
    "subgenus":         700,
    "species group":    600,
    "species subgroup": 500,
    "species":          400,
    "subspecies":       300,
    "infraspecies":     200,
    "varietas":         100,
    "forma":              0,
}
MINIMUM_SORTING_NUMBER = RANK_TO_SORTING_NUMBER['forma']
GENUS_RANK_TO_SORTING_NUMBER = RANK_TO_SORTING_NUMBER['genus']
SPECIES_SORTING_NUMBER = RANK_TO_SORTING_NUMBER['species']
MINIMUM_HIGHER_TAXON_NUMBER = 1 + SPECIES_SORTING_NUMBER
MAX_INFRASPECIFIC_NUMBER = SPECIES_SORTING_NUMBER - 1


def is_higher_taxon(r):
    return RANK_TO_SORTING_NUMBER[r] > SPECIES_SORTING_NUMBER


def is_infraspecific(r):
    return RANK_TO_SORTING_NUMBER[r] < SPECIES_SORTING_NUMBER


def is_species(r):
    return r == 'species'

