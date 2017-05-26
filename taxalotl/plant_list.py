#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR
from __future__ import print_function
import os

from peyotl import (get_logger, download_large_file)

from taxalotl.ott_schema import InterimTaxonomyData
from taxalotl.resource_wrapper import TaxonomyWrapper

_LOG = get_logger(__name__)


class PlantListWrapper(TaxonomyWrapper):
    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    def download(self):
        dd = self.unpacked_filepath
        _LOG.info("uf = {}".format(dd))
        top_files = []
        for u in self.url_list:
            pref, suff = os.path.split(u)
            if not suff:
                pref, suff = os.path.split(pref)
            _LOG.info("p = {} s = {}".format(pref, suff))
            assert suff
            dfp = os.path.join(dd, suff)
            top_files.append(dfp)
            if not os.path.exists(dfp):
                _LOG.debug("Starting download from {} to {}".format(u, dfp))
                download_large_file(u, dfp)
                _LOG.debug("Download from {} to {} completed.".format(u, dfp))
        for dfp in top_files:
            scrape_families_from_higher_group(dd, dfp)

    @property
    def download_filepath(self):
        return os.path.join(self.unpacked_filepath, "download_complete.txt")
