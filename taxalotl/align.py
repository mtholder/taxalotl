#!/usr/bin/env python
from __future__ import print_function

import sys
import os
from enum import IntEnum, IntFlag
from typing import Dict, List

from peyotl import (get_logger, write_as_json)

from .config import TaxalotlConfig
from .tree import TaxonTree
from .partitions import (PART_NAMES)
from .resource_wrapper import TaxonomyWrapper
from .taxonomic_ranks import (GENUS_RANK_TO_SORTING_NUMBER,
                              MAX_INFRASPECIFIC_NUMBER,
                              MINIMUM_HIGHER_TAXON_NUMBER,
                              SPECIES_SORTING_NUMBER)
from .tax_partition import IGNORE_SYN_TYPES

_LOG = get_logger(__name__)
out_stream = sys.stdout

def align_resource(taxalotl_config: TaxalotlConfig,
                   ott_res: TaxonomyWrapper,
                   res: TaxonomyWrapper,
                   level_list: List[str]):
    m = 'Could not align taxonomy because {} has not been partitioned.'
    for el in [ott_res, res]:
        if not el.has_been_partitioned():
            raise RuntimeError(m.format(el.id))
    if level_list == [None]:
        level_list = PART_NAMES
    for part_name in level_list:
        align_for_level(taxalotl_config, ott_res, res, part_name)

def align_for_level(taxalotl_config: TaxalotlConfig,
                    ott_res: TaxonomyWrapper,
                             res: TaxonomyWrapper,
                             part_name: str):
    fragment = taxalotl_config.get_fragment_from_part_name(part_name)
    print(fragment)
