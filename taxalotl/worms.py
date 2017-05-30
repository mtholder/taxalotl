#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR
from __future__ import print_function
from peyotl import (get_logger)
from taxalotl.ott_schema import InterimTaxonomyData
from taxalotl.resource_wrapper import TaxonomyWrapper

_LOG = get_logger(__name__)



class WormsTaxonomyWrapper(TaxonomyWrapper):
    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)
