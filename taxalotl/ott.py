from __future__ import print_function

import codecs
import os
from peyotl import (add_or_append_to_dict, get_logger, assure_dir_exists)
from taxalotl.partitions import INP_TAXONOMY_DIRNAME, MISC_DIRNAME
_LOG = get_logger(__name__)

OTT_PARTMAP = {'Archaea': frozenset([996421]),
               'Bacteria': frozenset([844192]),
               'Eukaryota/SAR': frozenset([5246039]),
               'Eukaryota/Haptophyta': frozenset([151014]),
               'Eukaryota/Rhodophyta': frozenset([878953]),
               'Eukaryota/Archaeplastida': frozenset([5268475]),
               'Eukaryota/Glaucophyta': frozenset([664970]),
               'Eukaryota/Chloroplastida': frozenset([361838]),
               'Eukaryota/Fungi': frozenset([352914]),
               'Eukaryota/Metazoa': frozenset([691846]),
               'Eukaryota/Metazoa/Annelida': frozenset([941620]),
               'Eukaryota/Metazoa/Arthropoda': frozenset([632179]),
               'Eukaryota/Metazoa/Arthropoda/Malacostraca': frozenset([212701]),
               'Eukaryota/Metazoa/Arthropoda/Arachnida': frozenset([511967]),
               'Eukaryota/Metazoa/Arthropoda/Insecta': frozenset([1062253]),
               'Eukaryota/Metazoa/Arthropoda/Insecta/Diptera': frozenset([661378]),
               'Eukaryota/Metazoa/Arthropoda/Insecta/Coleoptera': frozenset([865243]),
               'Eukaryota/Metazoa/Arthropoda/Insecta/Lepidoptera': frozenset([965954]),
               'Eukaryota/Metazoa/Arthropoda/Insecta/Hymenoptera': frozenset([753726]),
               'Eukaryota/Metazoa/Bryozoa': frozenset([442934]),
               'Eukaryota/Metazoa/Chordata': frozenset([125642]),
               'Eukaryota/Metazoa/Cnidaria': frozenset([641033]),
               'Eukaryota/Metazoa/Ctenophora': frozenset([641212]),
               'Eukaryota/Metazoa/Mollusca': frozenset([802117]),
               'Eukaryota/Metazoa/Nematoda': frozenset([395057]),
               'Eukaryota/Metazoa/Platyhelminthes': frozenset([555379]),
               'Eukaryota/Metazoa/Porifera': frozenset([67819]),
               'Viruses': frozenset([4807313]),
               }
OTT_3_SEPARATION_TAXA = OTT_PARTMAP
#cellular organisms	93302
# Eukaryota	304358

