#!/usr/bin/env python
from __future__ import print_function

import os
from peyotl import get_logger
from taxalotl.resource_manager import ResourceManager

_LOG = get_logger(__name__)


def _none_for_missing_config_get(config, section, option, default=None):
    try:
        return config.get(section, option)
    except Exception:
        return default


class TaxalotlConfig(object):
    def __init__(self,
                 filepath=None,
                 raw_downloads_dir=None,
                 processed_dir=None,
                 normalized_dir=None,
                 partitioned_dir=None):
        if filepath is None:
            if os.path.exists("taxalotl.conf"):
                filepath = os.path.abspath("taxalotl.conf")
            else:
                home_loc = os.path.expanduser("~/.opentreeoflife/taxalotl/taxalotl.conf")
                if os.path.exists(home_loc):
                    filepath = home_loc
        if filepath is None:
            m = "filepath to taxalotl.conf must be provided (or it must be in the current dir, " \
                "or ~/.opentreeoflife/taxalotl/taxalotl.conf)."
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
        resd = _none_for_missing_config_get(cfg, 'paths', 'resources')
        normd = normalized_dir
        if normd is None:
            normd = _none_for_missing_config_get(cfg, 'paths', 'normalized')
        partsd = partitioned_dir
        if partsd is None:
            partsd = _none_for_missing_config_get(cfg, 'paths', 'partitioned')
        self.raw_downloads_dir = rdd
        self.processed_dir = pd
        self.normalized_dir = normd
        self.partitioned_dir = partsd
        self.resources_dir = resd
        self._resources_mgr = None
        cws = _none_for_missing_config_get(cfg, 'behavior', 'crash_with_stacktraces')
        if cws:
            cws = cfg.getboolean('behavior', 'crash_with_stacktraces')
        self.crash_with_stacktraces = bool(cws)
        assert self.resources_mgr is not None

    @property
    def resources_mgr(self):
        if self._resources_mgr is None:
            if self.resources_dir is None:
                m = "The resources dir must be provided in the config file or the current dir"
                raise RuntimeError(m)
            if not os.path.isdir(self.resources_dir):
                m = 'The resources dir "{}" does not exist. Creating an empty one...'
                _LOG.warn(m.format(self.resources_dir))
                os.makedirs(self.resources_dir)
            self._resources_mgr = ResourceManager(self.resources_dir)
            for v in self._resources_mgr.resources.values():
                v.config = self
        return self._resources_mgr

    def get_resource_by_id(self, res_id):
        # from taxalotl.resource_wrapper import ResourceWrapper type: (str) -> ResourceWrapper
        try:
            return self.resources_mgr.resources[res_id]
        except Exception:
            raise ValueError("Unknown resource ID '{}'".format(res_id))

    def get_terminalized_res_by_id(self, res_id, logging_action_str=None):
        return self.get_all_terminalized_res_by_id(res_id)[-1]

    def get_all_terminalized_res_by_id(self, res_id):
        orig_res = self.get_resource_by_id(res_id)
        if orig_res.id == res_id and not orig_res.is_abstract:
            return [orig_res]
        base_res = self.get_resource_by_id(orig_res.base_id)
        if base_res.is_abstract:
            return list(base_res.get_self_and_children())
        return [orig_res]

    def get_all_ancs(self, res):
        all_ancs = [res]
        cr = res
        par_id = getattr(cr, 'inherits_from', None)
        while par_id is not None:
            cr = self.get_resource_by_id(par_id)
            all_ancs.append(cr)
            par_id = getattr(cr, 'inherits_from', None)
        all_ancs.reverse()
        return all_ancs

    def get_known_license_info(self, res):
        """Returns lists of the unique elements in license_url and license_or_tou_info"""
        aa = self.get_all_ancs(res)
        aa.reverse()
        urls = []
        uset = set()
        terms_of_use = []
        touset = set()
        for r in aa:
            u, t = r.license_url, r.license_or_tou_info
            if u and (u not in uset):
                urls.append(u)
                uset.add(u)
            if t and (t not in touset):
                touset.add(t)
                terms_of_use.append(t)
        return urls, terms_of_use
