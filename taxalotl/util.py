#!/usr/bin/env python
# from __future__ import print_function

import os
import io
from peyotl import get_logger

_LOG = get_logger(__name__)

INTERACTIVE_MODE = False


def get_true_false_repsonse(p, true_func=None, def_value=False):
    if not INTERACTIVE_MODE:
        _LOG.warn('non-interactive mode. Answering {} to "{}"'.format(def_value, p))
        return def_value
    if true_func is None:
        # noinspection PyPep8
        true_func = lambda r: r.lower() == 'y'
    try:
        resp = input(p)
        return true_func(resp)
    except:
        return def_value


def unlink(fp):
    _LOG.info('Removing "{}" ...'.format(fp))
    os.unlink(fp)


_FILES_WRITTEN = []


def get_filepaths_overwritten():
    return list(_FILES_WRITTEN)


class OutFile(object):
    def __init__(self, filepath, mode='w', encoding='utf-8'):
        self.filepath = filepath
        self.mode = mode
        self.encoding = encoding
        self.out_stream = None

    def __enter__(self):
        self.out_stream = io.open(self.filepath, mode=self.mode, encoding=self.encoding)
        _FILES_WRITTEN.append(self.filepath)
        return self.out_stream

    def __exit__(self):
        if self.out_stream is not None:
            self.out_stream.close()
            self.out_stream = None
