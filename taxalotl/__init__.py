#!/usr/bin/env python
from __future__ import print_function
from peyotl import get_logger
import codecs
import json
import os
_LOG = get_logger(__name__)

_ALLOWED_RESOURCE_KEYS = frozenset(["resources", "references", "ott_versions", "next"])
def read_resource_file(fp):
    try:
        with codecs.open(fp, 'rU', encoding='utf-8') as inp:
            o = json.load(inp)
        for k in o.keys():
            if k not in _ALLOWED_RESOURCE_KEYS:
                raise RuntimeError("Unrecognized key '{}'".format(k))
        ref = o.get("resources")
        if ref and isinstance(ref, list):
            rd = {}
            for el in ref:
                key = el["id"]
                if key in rd:
                    raise RuntimeError('The id "{}" was repeated'.format(key))
                rd[key] = el
            o["resources"] = rd
        return o
    except:
        _LOG.exception("Error reading JSON from \"{}\"".format(fp))
        raise

def write_resources_file(obj, fp):
    with codecs.open(fp, 'w', encoding='utf-8') as outp:
        json.dump(obj, outp, indent=2)

class ResourceManager(object):
    _MERGED_FILENAME = ".merged.json"
    def __init__(self, resource_dir, update_merged=True):
        self.resource_dir = resource_dir
        self.resources = {}
        self.references = {}
        self._filepath_read = None
        if update_merged:
            self._update_merged()
        self.read_merged()

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
                print(f)
                ni = os.path.join(self.resource_dir, f)
                inputs.append(ni)
                if (not needs_update):
                    if os.path.getmtime(ni) > mtt:
                        needs_update = True
        if needs_update:
            md = {}
            for fp in inputs:
                o = read_resource_file(fp)
                for r in ["references", "resources"]:
                    mrd = md.get(r, {})
                    mrk = set(mrd.keys())
                    ord = o.get(r, {})
                    ork = set(ord.keys())
                    i = ork.intersection(mrk)
                    if i:
                        m = 'IDs repeated in {} and previous resource files: "{}"'
                        raise RuntimeError(m.format(fp, '" "'.join(list(i))))
                    mrd.update(ord)
                    md[r] = mrd
                    print(mrd)
                for u in ["ott_versions", "next"]:
                    if u in o:
                        if u in md:
                            m = 'Unique key "{}" repeated in {} and previous resource files'
                            raise RuntimeError(m.format(u, fp))
                        else:
                            md[u] = o[u]
            write_resources_file(md, mfp)

    def read_merged(self):
        mfp = os.path.join(self.resource_dir, ResourceManager._MERGED_FILENAME)
        o = read_resource_file(mfp)
        self.resources.update(o.get("resources", {}))
        self.references.update(o.get("references", {}))
        print("resources = {}".format(self.resources))
        print("references = {}".format(self.references))
        self._filepath_read = mfp
