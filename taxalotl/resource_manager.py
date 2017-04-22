#!/usr/bin/env python
from __future__ import print_function
from peyotl import (assure_dir_exists,
                    download_large_file,
                    get_logger, gunzip, gunzip_and_untar,
                    write_as_json)
import codecs
import json
import time
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
        json.dump(obj, outp, indent=2, sort_keys=True, separators=(',', ': '))


_known_res_attr = frozenset(['aliases',
                             'copy_status',
                             'date',
                             'format',
                             'id',
                             'inherits_from',
                             'inputs',
                             'latest_download_url',
                             'maintainer',
                             'preceded_by',
                             'references',
                             'resource_type',
                             'schema',
                             'source',
                             'url',
                             'version'
                             ])


class _ResWrapper(object):
    def __init__(self, obj, parent=None, refs=None):
        for k in _known_res_attr:
            self.__dict__[k] = obj.get(k)
        for k in obj.keys():
            if k not in _known_res_attr:
                raise RuntimeError("Key '{}' not recognized".format(k))
        self.parent = parent
        self.children = []
        if self.parent:
            parent.children.append(self)
            for k in _known_res_attr:
                pv = getattr(parent, k, None)
                if obj.get(k) is None and pv is not None:
                    setattr(self, k, pv)
        self.id = obj['id']
        if self.references:
            x = [self.references] if not isinstance(self.references, list) else self.references
            for r in x:
                if r not in refs:
                    _LOG.warn(
                        "reference '{}' in resource '{}' was not recognized.".format(r, self.id))
        self.config = None

    def get_leaf_obj(self):
        if self.children:
            return self.children[-1].get_leaf_obj()
        return self

    def download_filepath(self, config=None):
        if self.is_abstract:
            return None
        fn = os.path.split(self.url)[-1]
        c = config if config else self.config
        return os.path.join(c.raw_downloads_dir, fn)

    def unpacked_filepath(self, config=None):
        if self.is_abstract:
            return None
        c = config if config else self.config
        return os.path.join(c.raw_downloads_dir, self.id)

    def normalized_filepath(self, config=None):
        if self.is_abstract:
            return None
        c = config if config else self.config
        return os.path.join(c.normalized_dir, self.id)

    @property
    def is_abstract(self):
        return self.format is None or self.url is None or self.schema is None

    def has_been_downloaded(self, config=None):
        dfp = self.download_filepath(config)
        return dfp is not None and os.path.exists(dfp)

    def has_been_unpacked(self, config=None):
        dfp = self.unpacked_filepath(config)
        return dfp is not None and os.path.exists(dfp)

    def has_been_normalized(self, config=None):
        dfp = self.normalized_filepath(config)
        return dfp is not None \
               and os.path.exists(dfp) \
               and os.path.exists(os.path.join(dfp, 'taxonomy.tsv'))

    def download(self, config=None):
        dfp = self.download_filepath(config)
        _LOG.debug("Starting download from {} to {}".format(self.url, dfp))
        download_large_file(self.url, dfp)
        _LOG.debug("Download from {} to {} completed.".format(self.url, dfp))

    def unpack(self, config):
        if self.format == 'tar+gzip':
            ufp = self.unpacked_filepath(config)
            dfp = self.download_filepath(config)
            _LOG.debug("Starting gunzip_and_untar from {} to {}".format(dfp, ufp))
            gunzip_and_untar(dfp, ufp)
            _LOG.debug("gunzip_and_untar from {} to {} completed.".format(dfp, ufp))
        else:
            raise NotImplementedError(
                "Unpacking from {} format is not currently supported".format(self.format))

    def normalize(self, config):
        schema = self.schema.lower()
        ufp = self.unpacked_filepath(config)
        nfp = self.normalized_filepath(config)
        if schema == "ncbi taxonomy":
            normalize_ncbi(ufp, nfp, self.url)
        else:
            raise NotImplementedError(
                "Unpacking from {} format is not currently supported".format(self.format))

    def write_status(self, out, config, indent=''):
        dfp = self.download_filepath(config)
        if dfp is None:
            lo = self.get_leaf_obj()
            if lo == self:
                suff = ''
            else:
                suff = ' the most recent version appears to be "{}"'.format(lo.id)
            template = "{}: {}. {}. This is a class of resource, not a downloadable artifact{}.\n"
            out.write(template.format(self.id, self.resource_type, self.source, suff))
            return
        out.write('{}:\n{}{}: {}'.format(self.id, indent, self.resource_type, self.source))
        if self.version:
            out.write(' version {}\n'.format(self.version))
        else:
            out.write(' (unversioned)\n')
        out.write("{}date: {}\n".format(indent, self.date if self.date else 'unknown'))
        s = "is at" if os.path.exists(dfp) else "not yet downloaded to"
        out.write("{}Raw ({} format) {} {}\n".format(indent, self.format, s, dfp))
        ufp = self.unpacked_filepath(config)
        s = "is at" if os.path.exists(ufp) else "not yet unpacked to"
        out.write("{}Raw ({} schema) {} {}\n".format(indent, self.schema, s, ufp))


class ExternalTaxonomyWrapper(_ResWrapper):
    resource_type = 'external taxonomy'

    def __init__(self, obj, parent=None, refs=None):
        _ResWrapper.__init__(self, obj, parent=parent, refs=refs)
        # print("ET obj = {}".format(obj))


class OTTaxonomyWrapper(_ResWrapper):
    resource_type = 'open tree taxonomy'

    def __init__(self, obj, parent=None, refs=None):
        _ResWrapper.__init__(self, obj, parent=parent, refs=refs)


class OTTaxonomyIdListWrapper(_ResWrapper):
    resource_type = 'open tree taxonomy idlist'

    def __init__(self, obj, parent=None, refs=None):
        _ResWrapper.__init__(self, obj, parent=parent, refs=refs)


_wrapper_types = [OTTaxonomyWrapper, ExternalTaxonomyWrapper, OTTaxonomyIdListWrapper, ]


def get_resource_wrapper(raw, refs, parent=None):
    rt = raw["resource_type"].lower()
    for wt in _wrapper_types:
        if rt == wt.resource_type:
            return wt(raw, parent=parent, refs=refs)
    raise RuntimeError("resource_type '{}' not recognized".format(rt))


def get_subclass_resource_wrapper(raw, known_dict, refs):
    par_key = raw["inherits_from"]
    par = known_dict[par_key]
    raw["resource_type"] = par.resource_type
    return get_resource_wrapper(raw, refs, parent=par)


def _wrap_resources(res, refs):
    rd = {}
    by_par = {}
    for k, v in res.items():
        if k != v["id"]:
            raise RuntimeError("Key and id field do not match ({} != {})".format(k, v["id"]))
        par = v.get('inherits_from')
        if par:
            by_par.setdefault(par, []).append((k, v))
        else:
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
    return rd, refs


class ResourceManager(object):
    _MERGED_FILENAME = ".merged.json"

    def __init__(self, resource_dir, update_merged=True):
        self.resource_dir = resource_dir
        self.resources = None
        self.references = None
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
                o = read_resource_file(fp)
                for r in ["references", "resources"]:
                    mrd = md.get(r, {})
                    mrk = set(mrd.keys())
                    ordict = o.get(r, {})
                    ork = set(ordict.keys())
                    i = ork.intersection(mrk)
                    if i:
                        m = 'IDs repeated in {} and previous resource files: "{}"'
                        raise RuntimeError(m.format(fp, '" "'.join(list(i))))
                    mrd.update(ordict)
                    md[r] = mrd
                for u in ["ott_versions", "next"]:
                    if u in o:
                        if u in md:
                            m = 'Unique key "{}" repeated in {} and previous resource files'
                            raise RuntimeError(m.format(u, fp))
                        else:
                            md[u] = o[u]
            write_resources_file(md, mfp)

    def _read_merged(self):
        mfp = os.path.join(self.resource_dir, ResourceManager._MERGED_FILENAME)
        o = read_resource_file(mfp)
        self.resources, self.references = _wrap_resources(o.get("resources", {}),
                                                          o.get("references", {}))
        self._filepath_read = mfp
