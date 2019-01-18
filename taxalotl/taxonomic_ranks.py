#!/usr/bin/env python
# -*- coding: utf-8 -*-

_RANK_TO_SORTING_NUMBER = {
    "superkingdom": 280,
    "kingdom": 270,
    "subkingdom": 260,
    "superphylum": 250,
    "phylum": 240,
    "subphylum": 230,
    "superclass": 220,
    "class": 210,
    "subclass": 200,
    "infraclass": 190,
    "cohort": 180,
    "superorder": 170,
    "order": 160,
    "suborder": 150,
    "infraorder": 140,
    "parvorder": 130,
    "superfamily": 120,
    "family": 110,
    "subfamily": 100,
    "tribe": 90,
    "subtribe": 80,
    "genus": 70,
    "subgenus": 60,
    "species group": 50,
    "species subgroup": 40,
    "species": 30,
    "subspecies": 20,
    "infraspecies": 15,
    "varietas": 10,
    "forma": 0,
}
MINIMUM_SORTING_NUMBER = _RANK_TO_SORTING_NUMBER['forma']
GENUS_RANK_TO_SORTING_NUMBER = _RANK_TO_SORTING_NUMBER['genus']
SPECIES_SORTING_NUMBER = _RANK_TO_SORTING_NUMBER['species']
MINIMUM_HIGHER_TAXON_NUMBER = 1 + SPECIES_SORTING_NUMBER
MAX_INFRASPECIFIC_NUMBER = SPECIES_SORTING_NUMBER - 1
