#!/usr/bin/env python
# from __future__ import print_function

import os
from peyotl import get_logger

_LOG = get_logger(__name__)

def get_true_false_repsonse(p, true_func=None, def_value=False):
    if true_func is None:
        true_func = lambda r: r.lower() == 'y'
    try:
        resp = input(p)
        return true_func(resp)
    except:
        return def_value


def unlink(fp):
    _LOG.info('Removing "{}" ...'.format(fp))
    os.unlink(fp)
