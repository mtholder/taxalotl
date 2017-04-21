#!/usr/bin/env python
from __future__ import print_function
from peyotl import get_logger
from taxalotl.resource_manager import ResourceManager

import os

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
            raise ValueError("filepath to taxalotl.conf must be provided (or it must be in the current dir, or ~/.taxalotl).")
        if not os.path.isfile(filepath):
            raise ValueError('No config file found at "{}"'.format(filepath))
        self._filepath = filepath
        try:
            # noinspection PyCompatibility
            from ConfigParser import SafeConfigParser
        except ImportError:
            # noinspection PyCompatibility,PyUnresolvedReferences
            from configparser import ConfigParser as SafeConfigParser  # pylint: disable=F0401
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
        self.raw_downloads_dir = rdd
        self.processed_dir = pd
        self.resources_dir = resd
        self._resources_mgr = None

    @property
    def resources_mgr(self):
        if self._resources_mgr is None:
            if self.resources_dir is None:
                raise RuntimeError("The resources dir must be provided in the config file or the current dir")
            if not os.path.isdir(self.resources_dir):
                raise RuntimeError('"{}" is not an existing resources dir.'.format(self.resources_dir))
            self._resources_mgr = ResourceManager(self.resources_dir)
        return self._resources_mgr
