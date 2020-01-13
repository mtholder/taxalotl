#!/usr/bin/env python
# from __future__ import print_function

import json
import os
import io
from peyotl import get_logger

_LOG = get_logger(__name__)

INTERACTIVE_MODE = True


def _startswith_y(r):
    return r.lower() == 'y'


def get_true_false_repsonse(p, true_func=_startswith_y, def_value=False):
    if not INTERACTIVE_MODE:
        _LOG.warn('non-interactive mode. Answering {} to "{}"'.format(def_value, p))
        return def_value
    try:
        resp = input(p)
        return true_func(resp)
    except:
        return def_value


def unlink(fp):
    _LOG.info('Removing "{}" ...'.format(fp))
    os.unlink(fp)


_FILES_WRITTEN = []


class TaxalotlHistory(object):
    def __init__(self):
        self.hist_filepath = os.path.expanduser("~/.taxalotl_history")
        self.hist_content = None

    def _read_hist(self):
        if self.hist_content is not None:
            return
        if os.path.isfile(self.hist_filepath):
            with io.open(self.hist_filepath, 'r', encoding='utf-8') as hout:
                self.hist_content = json.load(hout)

    def _write_hist(self):
        if not self.hist_content:
            return
        with io.open(self.hist_filepath, mode='w', encoding='utf-8') as outp:
            json.dump(self.hist_content, outp, indent=2)

    def add_virtual_command(self, name, res_id=None, level=None, wrote_files=None):
        if not wrote_files:
            return
        try:
            if self.hist_content is None:
                self._read_hist()
        except:
            _LOG.exception('Exception reading taxalotl history')
            pass
        if not self.hist_content:
            self.hist_content = []
        record = {'command': name}
        if res_id:
            record['res_id'] = res_id
        if level:
            record['level'] = level
        if wrote_files:
            record['wrote_files'] = wrote_files
        self.hist_content.append(record)
        try:
            self._write_hist()
        except:
            _LOG.exception('Exception writing taxalotl history')
            pass


_HISTORY_WRAPPER = TaxalotlHistory()


def get_filepaths_overwritten():
    return list(_FILES_WRITTEN)


def clear_filepaths_overwritten():
    del _FILES_WRITTEN[:]


class VirtCommand(object):
    def __init__(self, name, res_id=None, level=None):
        self.name = name
        self.res_id = res_id
        self.level = level
        self.th = _HISTORY_WRAPPER

    def __enter__(self):
        clear_filepaths_overwritten()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        fpo = get_filepaths_overwritten()
        if not fpo:
            return
        self.th.add_virtual_command(self.name,
                                    res_id=self.res_id,
                                    level=self.level,
                                    wrote_files=fpo)
        clear_filepaths_overwritten()


class OutDir(object):
    def __init__(self, filepath, mode='w', encoding='utf-8'):
        self.filepath = filepath

    def __enter__(self):
        if not os.path.exists(self.filepath):
            _LOG.info('Creating directory {}'.format(self.filepath))
            os.makedirs(self.filepath)
            _FILES_WRITTEN.append(self.filepath)
        return self.filepath

    def __exit__(self, exc_type, exc_value, traceback):
        pass

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

    def __exit__(self, exc_type, exc_value, traceback):
        if self.out_stream is not None:
            self.out_stream.close()
            self.out_stream = None
