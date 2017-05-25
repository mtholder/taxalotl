#!/usr/bin/env python
from __future__ import print_function
import urllib
import shutil
import codecs
import json
import os

from peyotl import (assure_dir_exists,
                    download_large_file,
                    get_logger, gunzip, gunzip_and_untar,
                    unzip)
from taxalotl.newick import normalize_newick
from taxalotl.ncbi import normalize_ncbi
from taxalotl.irmng import normalize_irmng
from taxalotl.silva import normalize_silva_taxonomy, partition_silva
from taxalotl.darwin_core import normalize_darwin_core_taxonomy
from taxalotl.col import partition_col, partition_col_by_root_id
from taxalotl.ott import (ott_diagnose_new_separators,
                          ott_build_paritition_maps,
                          partition_ott_by_root_id,
                          )
from taxalotl.partitions import (find_partition_dirs_for_taxonomy,
                                 get_part_inp_taxdir,
                                 get_par_and_par_misc_taxdir,
                                 get_inp_taxdir,
                                 get_misc_inp_taxdir)
from taxalotl.interim_taxonomy_struct import (INP_OTT_SYNONYMS_HEADER,
                                              INP_OTT_TAXONOMY_HEADER)

_LOG = get_logger(__name__)


def get_resource_wrapper(raw, refs, parent=None):
    from taxalotl.resource_mapper import BASE_ID_TO_RES_TYPE, wrapper_types
    base = raw["base_id"].lower()
    if base:
        base(raw, parent=parent, refs)
    else:
        rt = raw["resource_type"].lower()
        for wt in wrapper_types:
            if rt == wt.resource_type:
                return wt(raw, parent=parent, refs=refs)
    raise RuntimeError("resource_type '{}' not recognized".format(rt))


def get_subclass_resource_wrapper(raw, known_dict, refs):
    par_key = raw["inherits_from"]
    par = known_dict[par_key]
    raw["resource_type"] = par.resource_type
    raw["base_id"] = par.base_id
    return get_resource_wrapper(raw, refs, parent=par)


def wrap_otifact_resources(res, refs=None):
    rd = {}
    by_par = {}
    for k, v in res.items():
        v["id"] = k
        par = v.get('inherits_from')
        if par:
            by_par.setdefault(par, []).append((k, v))
        else:
            v["base_id"] = k
            w = get_resource_wrapper(v, refs)
            rd[k] = w
    while by_par:
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
                rd[rk] = get_subclass_resource_wrapper(v, rd, refs)
    aliases = {}
    for w in rd.values():
        if w.aliases:
            for a in w.aliases:
                if a in aliases or a in rd:
                    raise RuntimeError("Repeated alias for an id: {}".format(a))
                aliases[a] = w
    rd.update(aliases)
    res_list = rd.values()
    for res in res_list:
        if res.base_id == 'col':
            res.__class__ = CoLExternalTaxonomyWrapper
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
