#!/usr/bin/env python
from __future__ import print_function

import os
from peyotl import get_logger
from taxalotl.resource_manager import ResourceManager
from taxalotl.resource_manager import ResourceWrapper

_LOG = get_logger(__name__)


def _none_for_missing_config_get(config, section, option, default=None):
    try:
        return config.get(section, option)
    except:
        return default


class TaxalotlConfig(object):
    def __init__(self,
                 filepath=None,
                 raw_downloads_dir=None,
                 processed_dir=None,
                 normalized_dir=None,
                 resources_dir=None):
        def_resources = ''
        if filepath is None:
            if os.path.exists("taxalotl.conf"):
                filepath = os.path.abspath("taxalotl.conf")
                if resources_dir is None and os.path.exists('resources'):
                    def_resources = os.path.abspath("resources")
            else:
                home_loc = os.path.expanduser("~/.taxalotl")
                if os.path.exists(home_loc):
                    filepath = home_loc
        if filepath is None:
            m = "filepath to taxalotl.conf must be provided (or it must be in the current dir, " \
                "or ~/.taxalotl)."
            raise ValueError(m)
        if not os.path.isfile(filepath):
            raise ValueError('No config file found at "{}"'.format(filepath))
        self._filepath = filepath
        try:
            # noinspection PyCompatibility,PyUnresolvedReferences
            from ConfigParser import SafeConfigParser
        except ImportError:
            # noinspection PyCompatibility,PyUnresolvedReferences
            from configparser import ConfigParser as SafeConfigParser  # pylint: disable=F0401
        # noinspection PyUnboundLocalVariable
        cfg = SafeConfigParser()
        cfg.read([filepath])
        #
        rdd = raw_downloads_dir
        if rdd is None:
            rdd = _none_for_missing_config_get(cfg, 'paths', 'raw')
        pd = processed_dir
        if pd is None:
            pd = _none_for_missing_config_get(cfg, 'paths', 'processed')
        resd = resources_dir
        if resd is None:
            resd = _none_for_missing_config_get(cfg, 'paths', 'resources', def_resources)
        normd = normalized_dir
        if normd is None:
            normd = _none_for_missing_config_get(cfg, 'paths', 'normalized')
        self.raw_downloads_dir = rdd
        self.processed_dir = pd
        self.normalized_dir = normd
        self.resources_dir = resd
        self._resources_mgr = None
        cws = _none_for_missing_config_get(cfg, 'behavior', 'crash_with_stacktraces')
        if cws:
            cws = cfg.getboolean('behavior', 'crash_with_stacktraces')
        self.crash_with_stacktraces = bool(cws)

    @property
    def resources_mgr(self):
        if self._resources_mgr is None:
            if self.resources_dir is None:
                raise RuntimeError(
                    "The resources dir must be provided in the config file or the current dir")
            if not os.path.isdir(self.resources_dir):
                raise RuntimeError(
                    '"{}" is not an existing resources dir.'.format(self.resources_dir))
            self._resources_mgr = ResourceManager(self.resources_dir)
            for v in self._resources_mgr.resources.values():
                v.config = self
        return self._resources_mgr

    def get_resource_by_id(self, res_id):
        # type: (str) -> ResourceWrapper
        try:
            return self.resources_mgr.resources[res_id]
        except:
            raise ValueError("Unknown resource ID '{}'".format(res_id))

    def get_terminalized_res_by_id(self, res_id, logging_action_str=None):
        orig_res = self.get_resource_by_id(res_id)
        if orig_res.is_abstract:
            wrapper = orig_res.get_leaf_obj()
            if logging_action_str and (wrapper is not orig_res):
                m = "{} action being performed on {} as the most recent version of {}"
                _LOG.info(m.format(logging_action_str, wrapper.id, orig_res.id))
            return wrapper
        return orig_res
