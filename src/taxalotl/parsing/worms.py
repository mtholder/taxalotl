#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR
from __future__ import print_function
import logging


from ..resource_wrapper import TaxonomyWrapper

_LOG = logging.getLogger(__name__)


class WormsTaxonomyWrapper(TaxonomyWrapper):
    schema = {"ott"}

    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)
