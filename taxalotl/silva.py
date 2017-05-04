#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR and from:
#   reference-taxonomy/feed/gbif/project_2016.py
# and
#   reference-taxonomy/feed/gbif/process_gbif_taxonomy.py
from __future__ import print_function

import codecs
import os
import re

from peyotl import (assure_dir_exists,
                    get_logger)

from taxalotl.interim_taxonomy_struct import InterimTaxonomyData

_LOG = get_logger(__name__)

def normalize_silva_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    itd = InterimTaxonomyData()
    sys.exit('hi from normalize silva')