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
from taxalotl.darwin_core import normalize_darwin_core_taxonomy

_LOG = get_logger(__name__)

_ALLOWED_RESOURCE_KEYS = frozenset(["resources", "references", "ott_versions", "next"])


######################################################################################
def unpack_archive(archive_fp, unpack_fp, archive_format, wrapper):
    afl = archive_format.lower()
    if afl in ['tar+gzip']:
        _LOG.debug("gunzip_and_untar from {} to {} ...".format(archive_fp, unpack_fp))
        gunzip_and_untar(archive_fp, unpack_fp)
        _LOG.debug("gunzip_and_untar from {} to {} done.".format(archive_fp, unpack_fp))
    elif afl == 'zip':
        _LOG.debug("unzip from {} to {} ...".format(archive_fp, unpack_fp))
        unzip(archive_fp, unpack_fp)
        _LOG.debug("unzip from {} to {} done.".format(archive_fp, unpack_fp))
    elif afl == 'gzip':
        afn = os.path.split(archive_fp)[-1]
        if archive_fp.endswith(".gz"):
            fn = afn[:-3]
        elif archive_fp.endswith(".gzip"):
            fn = afn[:-5]
        else:
            raise RuntimeError("Expecting gzipped archive to endwith .gz or .gzip")
        assure_dir_exists(unpack_fp)
        dest = os.path.join(unpack_fp, fn)
        _LOG.debug("gunzip from {} to {} ...".format(archive_fp, dest))
        gunzip(archive_fp, dest)
        _LOG.debug("gunzip from {} to {} done.".format(archive_fp, dest))
    elif afl == 'text':
        assure_dir_exists(unpack_fp)
        try:
            lfn = wrapper.local_filename
            assert lfn is not None
        except:
            raise RuntimeError("Resource must have a local_filename if it format=text")
        shutil.copyfile(archive_fp, os.path.join(unpack_fp, lfn))
    else:
        m = "Unpacking from {} format is not currently supported"
        raise NotImplementedError(m.format(archive_format))


######################################################################################
OTT_TAXONOMY_FILENAMES = ("about.json",
                          "alignment_summary.json",
                          "choices.tsv",
                          "conflicts.tsv",
                          "deprecated.tsv",
                          "differences.tsv",
                          "forwards.tsv",
                          "log.tsv",
                          "merge_summary.json",
                          "otu_differences.tsv",
                          "README.html",
                          "synonyms.tsv",
                          "taxonomy.tsv",
                          "transcript.out",
                          "version.txt",)

OTT_TAXONOMY_ID_FILES = ('by_qid.csv', )

def copy_file_list_by_linking(unpacked_dirp, normalized_dirp, file_list):
    assure_dir_exists(normalized_dirp)
    for fn in file_list:
        ufp = os.path.join(unpacked_dirp, fn)
        if os.path.exists(ufp):
            dfp = os.path.join(normalized_dirp, fn)
            if os.path.exists(dfp):
                _LOG.info('File already exists at "{}". Skipping link creation.'.format(dfp))
            else:
                os.symlink(ufp, dfp)


def copy_taxonomy_by_linking(unpacked_dirp, normalized_dirp, resource_wrapper):
    copy_file_list_by_linking(unpacked_dirp,
                              normalized_dirp,
                              OTT_TAXONOMY_FILENAMES)

def copy_id_list_by_linking(unpacked_dirp, normalized_dirp, resource_wrapper):
    copy_file_list_by_linking(unpacked_dirp,
                              normalized_dirp,
                              OTT_TAXONOMY_ID_FILES)


_schema_to_norm_fn = {"ott": copy_taxonomy_by_linking,
                      "ott id csv": copy_id_list_by_linking,
                      "ncbi taxonomy": normalize_ncbi,
                      "http://rs.tdwg.org/dwc/": normalize_darwin_core_taxonomy,
                      "newick": normalize_newick,
                      "irmng dwc": normalize_irmng,
                      }


def normalize_archive(unpacked_fp, normalized_fp, schema_str, resource_wrapper):
    schema = schema_str.lower()
    try:
        norm_fn = _schema_to_norm_fn[schema]
    except KeyError:
        m = "Normalization from \"{}\" schema is not currently supported"
        raise NotImplementedError(m.format(schema_str))
    norm_fn(unpacked_fp, normalized_fp, resource_wrapper)


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
                             'local_filename',
                             'maintainer',
                             'preceded_by',
                             'references',
                             'resource_type',
                             'schema',
                             'source',
                             'url',
                             'version'
                             ])


class ResourceWrapper(object):
    def __init__(self, obj, parent=None, refs=None, config=None):
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
        self._config = config

    @property
    def config(self):
        if self._config is None:
            m = "Inadequately initialized ResourceWrapper: " \
                "config attribute of a resource must be set."
            raise RuntimeError(m)
        return self._config

    @config.setter
    def config(self, c):
        self._config = c

    def get_leaf_obj(self):
        # type: () -> ResourceWrapper
        if self.children:
            return self.children[-1].get_leaf_obj()
        return self

    @property
    def download_filepath(self):
        if self.is_abstract:
            return None
        fn = os.path.split(self.url)[-1]
        return os.path.join(self.config.raw_downloads_dir, fn)

    @property
    def unpacked_filepath(self):
        if self.is_abstract:
            return None
        return os.path.join(self.config.raw_downloads_dir, self.id)

    @property
    def normalized_filepath(self):
        if self.is_abstract:
            return None
        return os.path.join(self.config.normalized_dir, self.id)

    @property
    def is_abstract(self):
        return self.format is None or self.url is None or self.schema is None

    def has_been_downloaded(self):
        dfp = self.download_filepath
        return dfp is not None and os.path.exists(dfp)

    def has_been_unpacked(self):
        dfp = self.unpacked_filepath
        return dfp is not None and os.path.exists(dfp)

    def has_been_normalized(self):
        dfp = self.normalized_filepath
        return (dfp is not None
                and os.path.exists(dfp)
                and os.path.exists(os.path.join(dfp, 'taxonomy.tsv')))

    def download(self):
        dfp = self.download_filepath
        if dfp is None:
            m = "Resource {} appears to be abstract, therefore not downloadable"
            raise RuntimeError(m.format(self.id))
        if self.url.startswith("ftp://"):
            _LOG.debug("Starting FTP download from {} to {}".format(self.url, dfp))
            urllib.urlretrieve(self.url, dfp)
            _LOG.debug("Download from {} to {} completed.".format(self.url, dfp))
        else:
            _LOG.debug("Starting download from {} to {}".format(self.url, dfp))
            download_large_file(self.url, dfp)
            _LOG.debug("Download from {} to {} completed.".format(self.url, dfp))

    def unpack(self):
        unpack_archive(self.download_filepath, self.unpacked_filepath, self.format, self)

    def normalize(self):
        normalize_archive(self.unpacked_filepath, self.normalized_filepath, self.schema, self)

    def write_status(self, out, indent='', list_all_artifacts=False):
        dfp = self.download_filepath
        if dfp is None:
            lo = self.get_leaf_obj()
            if lo == self:
                suff = ''
            else:
                suff = ' the most recent version appears to be "{}"'.format(lo.id)
            template = "{}: {}. {}. This is a class of resource, not a downloadable artifact{}.\n"
            out.write(template.format(self.id, self.resource_type, self.source, suff))
            return
        out.write('{}: {} {} '.format(self.id, self.resource_type, self.source))
        if self.version:
            out.write('version {}. '.format(self.version))
        else:
            out.write('(unversioned). ')
        out.write("date={}\n".format(self.date if self.date else 'unknown'))

        s = "is at" if self.has_been_downloaded() else "not yet downloaded to"
        down_str = "{}Raw ({} format) {} {}\n".format(indent, self.format, s, dfp)
        ufp = self.unpacked_filepath
        s = "is at" if self.has_been_unpacked() else "not yet unpacked to"
        unp_str = "{}Raw ({} schema) {} {}\n".format(indent, self.schema, s, ufp)
        nfp = self.normalized_filepath
        s = "is at" if self.has_been_normalized() else "not yet normalized to"
        norm_str = "{}OTT formatted form {} {}\n".format(indent, s, nfp)
        if list_all_artifacts:
            out.write('{}{}{}'.format(down_str, unp_str, norm_str))
        else:
            if self.has_been_normalized():
                out.write(norm_str)
            elif self.has_been_unpacked():
                out.write(unp_str)
            else:
                out.write(down_str)


# noinspection PyAbstractClass
class ExternalTaxonomyWrapper(ResourceWrapper):
    resource_type = 'external taxonomy'

    def __init__(self, obj, parent=None, refs=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs)
        # print("ET obj = {}".format(obj))


# noinspection PyAbstractClass
class OTTaxonomyWrapper(ResourceWrapper):
    resource_type = 'open tree taxonomy'

    def __init__(self, obj, parent=None, refs=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs)


# noinspection PyAbstractClass
class OTTaxonomyIdListWrapper(ResourceWrapper):
    resource_type = 'open tree taxonomy idlist'

    def __init__(self, obj, parent=None, refs=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs)

    def has_been_normalized(self):
        dfp = self.normalized_filepath
        return (dfp is not None
                and os.path.exists(dfp)
                and os.path.exists(os.path.join(dfp, 'by_qid.csv')))


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
