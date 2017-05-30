#!/usr/bin/env python
from __future__ import print_function
from peyotl import get_logger
import os
from taxalotl.tax_partition import INP_TAXONOMY_DIRNAME, MISC_DIRNAME

_LOG = get_logger(__name__)

def compare_taxonomies_in_dir(taxalotl_conf, tax_dir):
    tax_id_set = set()
    non_misc_dir = os.path.join(tax_dir, INP_TAXONOMY_DIRNAME)
    misc_dir = os.path.join(tax_dir, MISC_DIRNAME, INP_TAXONOMY_DIRNAME)
    for sd in [misc_dir, non_misc_dir]:
        if os.path.exists(sd):
            tax_id_set.update(os.listdir(sd))
    _LOG.info("tax_id = {}".format(tax_id_set))