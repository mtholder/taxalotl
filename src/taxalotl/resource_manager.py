#!/usr/bin/env python
from __future__ import print_function

import io
import json
import os
from .resource_wrapper import AbstractResourceWrapper
from .util import OutFile
import logging

_LOG = logging.getLogger(__name__)


def read_resource_file(fp):
    try:
        with io.open(fp, 'rU', encoding='utf-8') as inp:
            return json.load(inp)
    except:
        _LOG.exception("Error reading JSON from \"{}\"".format(fp))
        raise


def write_resources_file(obj, fp):
    with OutFile(fp) as outp:
        json.dump(obj, outp, indent=2, sort_keys=True, separators=(',', ': '))


def get_resource_wrapper(raw, refs, parent=None):
    from .resource_mapper import BASE_ID_TO_RES_TYPE, wrapper_types
    base = BASE_ID_TO_RES_TYPE.get(raw["base_id"].lower())
    if base and base is not AbstractResourceWrapper:
        # _LOG.debug('get_resource_wrapper.calling base raw={}'.format(raw))
        return base(raw, parent=parent, refs=refs)
    rt = raw["resource_type"].lower()
    st = raw.get('schema', '').lower()
    # _LOG.debug('rt={}'.format(rt))
    for wt in wrapper_types:
        if rt == wt.resource_type and (st in wt.schema):
            # _LOG.debug('get_resource_wrapper.calling wrapper_types wt={}'.format(wt))
            return wt(raw, parent=parent, refs=refs)
    # m = "resource_type, schema = ({}, {}) not recognized for {}, using AbstractResourceWrapper..." \
    #     .format(rt, st, raw['id'])
    # _LOG.info(m)
    return AbstractResourceWrapper(raw, parent=parent, refs=refs)


def get_subclass_resource_wrapper(raw, known_dict, refs):
    par_key = raw["inherits_from"]
    par = known_dict[par_key]
    raw["resource_type"] = par.resource_type
    raw["base_id"] = par.base_id
    # _LOG.debug('get_subclass_resource_wrapper.raw = {}'.format(raw))
    return get_resource_wrapper(raw, refs, parent=par)


def wrap_otifact_resources(res, refs=None):
    rd = {}
    by_par = {}
    for k, v in res.items():
        # _LOG.debug('k,v = {}, {}'.format(k, v))
        v["id"] = k
        par = v.get('inherits_from')
        if par:
            by_par.setdefault(par, []).append((k, v))
        else:
            v["base_id"] = k
            # _LOG.debug('wrap_otifact_resources k = {}'.format(k))
            w = get_resource_wrapper(v, refs)
            if w is not None:
                rd[k] = w

    while by_par:
        # for k, v in by_par.items():
        #     for rk, inner_v in v:
        #         _LOG.debug('Round {}: by_par[{}] = [{}, {}]'.format(round, k, rk, inner_v))
        # round += 1
        curr_key_set = set(rd.keys())
        needed = set(by_par.keys())
        inter = curr_key_set.intersection(needed)
        if not inter:
            raise RuntimeError(
                "Could not find the base class resources '{}'".format("', '").join(by_par.keys()))
        for k in inter:
            rk_v_list = by_par[k]
            del by_par[k]
            for rk, v in rk_v_list:
                # _LOG.debug('rk,v = {}, {}'.format(rk, v))
                rd[rk] = get_subclass_resource_wrapper(v, rd, refs)
    # _LOG.warning('rd = {}'.format(rd))
    aliases = {}
    for w in rd.values():
        if w.aliases:
            for a in w.aliases:
                pra = aliases.get(a)
                if pra:
                    if pra.base_id != pra.id:
                        aliases[a] = w
                if a in rd:
                    raise RuntimeError("Previously registered for an id: {}".format(a))
                aliases[a] = w
    # _LOG.warning('aliases = {}'.format(aliases))
    rd.update(aliases)
    return rd


class ResourceManager(object):
    _MERGED_FILENAME = ".merged.json"

    def __init__(self, resource_dir, update_merged=True):
        self.resource_dir = resource_dir
        self.resources = None
        self._filepath_read = None
        if update_merged:
            self._update_merged()
        self._read_merged()
        self._generic_handler = None

    def generic_handler_can_be_used(self, res_id):
        return res_id.startswith('cof-')

    def _update_merged(self):
        needs_update = False
        inputs = []
        mfp = os.path.join(self.resource_dir, ResourceManager._MERGED_FILENAME)
        if os.path.exists(mfp):
            mtt = os.path.getmtime(mfp)
        else:
            needs_update = True
        for f in os.listdir(self.resource_dir):
            if f.endswith(".json") and f != ResourceManager._MERGED_FILENAME:
                ni = os.path.join(self.resource_dir, f)
                inputs.append(ni)
                if not needs_update:
                    # noinspection PyUnboundLocalVariable
                    if os.path.getmtime(ni) > mtt:
                        needs_update = True
        if needs_update:
            md = {}
            for fp in inputs:
                mrk = set(md.keys())
                ordict = read_resource_file(fp)
                ork = set(ordict.keys())
                i = ork.intersection(mrk)
                if i:
                    m = 'IDs repeated in {} and previous resource files: "{}"'
                    raise RuntimeError(m.format(fp, '" "'.join(list(i))))
                md.update(ordict)
            write_resources_file(md, mfp)

    def _read_merged(self):
        mfp = os.path.join(self.resource_dir, ResourceManager._MERGED_FILENAME)
        self.resources = wrap_otifact_resources(read_resource_file(mfp))
        self._filepath_read = mfp

    def abstract_input_resource_types(self):
        return [k for k, v in self.resources.items() if v.is_abstract_input_resource_type]
