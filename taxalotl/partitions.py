#!/usr/bin/env python
# from __future__ import print_function


TOP_PARTS = ('Archaea',
             'Bacteria',
             'Eukaryota/Metazoa',
             'Eukaryota/Fungi',
             'Eukaryota/plants',
             'Eukaryota/Archaeplastida',
             'Eukaryota/__misc__',
             'Eukaryota/__other__',
             'Viruses',
             '__misc__',
             )

METAZOA_PARTS = ('Annelida'
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
METAZOA_FRAG = 'Eukaryota/Metazoa/'
METAZOA_PARTS  = tuple([METAZOA_FRAG + i for i in METAZOA_PARTS])
PARTS_BY_NAME = {'Life': TOP_PARTS,
                 'Metazoa': METAZOA_PARTS,
                }
PART_FRAG_BY_NAME= {'Life': '',
                    'Metazoa': METAZOA_FRAG,
                   }
PART_NAMES = list(PARTS_BY_NAME.keys())
PART_NAMES.sort()
PART_NAMES = tuple(PART_NAMES)