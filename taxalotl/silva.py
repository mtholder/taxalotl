#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Previous silva handling code that served as a basis for this code was written by JAR and
#   Jessica Grant as a part of the reference_taxonomy and OToL efforts.
from __future__ import print_function

from peyotl import (assure_dir_exists,
                    get_logger)
import codecs
import os
from taxalotl.interim_taxonomy_struct import InterimTaxonomyData
from taxalotl.commands import unpack_resources
_LOG = get_logger(__name__)

def parse_silva_ids(fn):
    id_map = {}
    with codecs.open(fn, 'r', encoding='utf-8') as inp:
        for line in inp:
            ls = line.strip()
            lsslash = ls.split('/')
            if ls in id_map:
                _LOG.warn("repeated ID: {} in {}".format(ls, fn))
            id_map[ls] = lsslash
    return id_map

def normalize_silva_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    itd = InterimTaxonomyData()
    id_list_id = res_wrapper.id_list
    taxalotl_config = res_wrapper.config
    id_list_res = taxalotl_config.get_terminalized_res_by_id(id_list_id, 'normalize silva')
    if not id_list_res.has_been_unpacked():
        unpack_resources(taxalotl_config, [id_list_id])
    expect_id_fp = os.path.join(id_list_res.unpacked_filepath, id_list_res.local_filename)
    if not os.path.isfile(expect_id_fp):
        raise ValueError("Silva ID file not found at: {}". format(expect_id_fp))
    preferred_ids = parse_silva_ids(expect_id_fp)