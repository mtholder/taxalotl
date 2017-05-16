#!/usr/bin/env python
# from __future__ import print_function

INP_TAXONOMY_DIRNAME = '__inputs__'
MISC_DIRNAME = '__misc__'
EUK_MICROBE_DIRNAME = '__other__'
TOP_PARTS = ('Archaea',
             'Bacteria',
             'Eukaryota/Metazoa',
             'Eukaryota/Fungi',
             'Eukaryota/plants',
             'Eukaryota/Archaeplastida',
             'Eukaryota/' + MISC_DIRNAME,
             'Eukaryota/' + EUK_MICROBE_DIRNAME,
             'Viruses',
             '__misc__',
             )

_part_list = [('Life', TOP_PARTS, ''),]
METAZOA_PARTS = ('Annelida',
                 'Arthropoda',
                 'Bryozoa',
                 'Chordata',
                 'Cnidaria',
                 'Ctenophora',
                 'Mollusca',
                 'Nematoda',
                 'Platyhelminthes',
                 'Porifera',
                 '__misc__',
                 )
METAZOA = 'Metazoa'
METAZOA_FRAG = 'Eukaryota/{}/'.format(METAZOA)
METAZOA_PARTS  = tuple([METAZOA_FRAG + i for i in METAZOA_PARTS])
_part_list.append((METAZOA, METAZOA_PARTS, METAZOA_FRAG))

ARTHROPODA = 'Arthropoda'
ARTHROPODA_PARTS = ('Malacostraca', 'Arachnida', 'Insecta')
ARTHROPODA_FRAG = METAZOA_FRAG + ARTHROPODA + '/'
ARTHROPODA_PARTS  = tuple([ARTHROPODA_FRAG + i for i in ARTHROPODA_PARTS])
_part_list.append((ARTHROPODA, ARTHROPODA_PARTS, ARTHROPODA_FRAG))

INSECTA = 'Insecta'
INSECTA_PARTS = ('Diptera', 'Coleoptera', 'Lepidoptera', 'Hymenoptera')
INSECTA_FRAG = ARTHROPODA_FRAG + INSECTA + '/'
INSECTA_PARTS  = tuple([INSECTA_FRAG + i for i in INSECTA_PARTS])
_part_list.append((INSECTA, INSECTA_PARTS, INSECTA_FRAG))

PARTS_BY_NAME = {}
PART_FRAG_BY_NAME = {}
for key, parts, frag in _part_list:
    PARTS_BY_NAME[key] = parts
    PART_FRAG_BY_NAME[key] = frag
PART_NAMES = list(PARTS_BY_NAME.keys())
PART_NAMES.sort()
PART_NAMES = tuple(PART_NAMES)