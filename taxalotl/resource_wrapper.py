#!/usr/bin/env python
from __future__ import print_function

import codecs
import os
import shutil
import urllib

from peyotl import (assure_dir_exists,
                    download_large_file,
                    get_logger, gunzip, gunzip_and_untar,
                    unzip)

from taxalotl.ott_schema import (INP_OTT_SYNONYMS_HEADER,
                                 INP_OTT_TAXONOMY_HEADER,
                                 partition_ott_by_root_id)
from taxalotl.newick import normalize_newick
from taxalotl.partitions import (find_partition_dirs_for_taxonomy,
                                 get_auto_gen_part_mapper,
                                 get_part_inp_taxdir,
                                 get_par_and_par_misc_taxdir,
                                 get_inp_taxdir,
                                 get_misc_inp_taxdir,
                                 do_partition)

_LOG = get_logger(__name__)


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
        if wrapper.local_filename:
            dest = os.path.join(unpack_fp, wrapper.local_filename)
        else:
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

OTT_TAXONOMY_ID_FILES = ('by_qid.csv',)


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


# noinspection PyUnusedLocal
def copy_taxonomy_by_linking(unpacked_dirp, normalized_dirp, resource_wrapper):
    copy_file_list_by_linking(unpacked_dirp,
                              normalized_dirp,
                              OTT_TAXONOMY_FILENAMES)


# noinspection PyUnusedLocal
def copy_id_list_by_linking(unpacked_dirp, normalized_dirp, resource_wrapper):
    copy_file_list_by_linking(unpacked_dirp,
                              normalized_dirp,
                              OTT_TAXONOMY_ID_FILES)


# noinspection PyUnusedLocal
def copy_and_add_ott_headers(unpacked_dirp, normalized_dirp, resource_wrapper):
    motf = list(OTT_TAXONOMY_FILENAMES)
    special = [('taxonomy.tsv', INP_OTT_TAXONOMY_HEADER),
               ('synonyms.tsv', INP_OTT_SYNONYMS_HEADER)]
    for fn, header in special:
        motf.remove(fn)
    copy_file_list_by_linking(unpacked_dirp, normalized_dirp, motf)
    for fn, header in special:
        tf = os.path.join(unpacked_dirp, fn)
        if os.path.isfile(tf):
            content = codecs.open(tf, 'r', encoding='utf-8').read()
            outfp = os.path.join(normalized_dirp, fn)
            with codecs.open(outfp, 'w', encoding='utf-8') as out:
                out.write(header)
                out.write(content)


# noinspection PyUnusedLocal
def normalize_tab_sep_ott(unpacked_dirp, normalized_dirp, resource_wrapper):
    motf = list(OTT_TAXONOMY_FILENAMES)
    special = [('taxonomy.tsv', INP_OTT_TAXONOMY_HEADER)]
    for fn, header in special:
        motf.remove(fn)
    copy_file_list_by_linking(unpacked_dirp, normalized_dirp, motf)
    for fn, header in special:
        tf = os.path.join(unpacked_dirp, fn)
        if os.path.isfile(tf):
            outfp = os.path.join(normalized_dirp, fn)
            with codecs.open(tf, 'r', encoding='utf-8') as inp:
                with codecs.open(outfp, 'w', encoding='utf-8') as out:
                    for line in inp:
                        ls = line.split('\t')
                        out.write('\t|\t'.join(ls))


_schema_to_norm_fn = {"headerless ott": copy_and_add_ott_headers,
                      "newick": normalize_newick,
                      "ott": copy_taxonomy_by_linking,
                      "ott id csv": copy_id_list_by_linking,
                      "tab-separated ott": normalize_tab_sep_ott,
                      }

_known_res_attr = frozenset(['aliases',
                             'base_id',  # base ID in inherits from graph
                             'copy_status',
                             'date',
                             'depends_on',
                             'doc_url',
                             'format',
                             'id',
                             'inherits_from',
                             'inputs',
                             'id_list',
                             'latest_download_url',
                             'license_url',
                             'license_or_tou_info',
                             'local_filename',
                             'maintainer',
                             'notes',
                             'preceded_by',
                             'references',
                             'resource_type',
                             'schema',
                             'source',
                             'stats',
                             'url',
                             'url_list',
                             'version'
                             ])


class FromOTifacts(object):
    def __init__(self):
        self.aliases = None
        self.base_id = None
        self.copy_status = None
        self.date = None
        self.depends_on = None
        self.doc_url = None
        self.format = None
        self.id = None
        self.inherits_from = None
        self.inputs = None
        self.id_list = None
        self.latest_download_url = None
        self.license_url = None
        self.license_or_tou_info = None
        self.local_filename = None
        self.maintainer = None
        self.notes = None
        self.preceded_by = None
        self.references = None
        self.resource_type = None
        self.schema = None
        self.source = None
        self.stats = None
        self.url = None
        self.url_list = None
        self.version = None


class ResourceWrapper(FromOTifacts):
    taxon_filename = 'taxonomy.tsv'
    synonyms_filename = 'synonyms.tsv'
    partition_parsing_fn = staticmethod(partition_ott_by_root_id)

    def __init__(self, obj, parent=None, refs=None, config=None):
        FromOTifacts.__init__(self)
        self.base_id = None
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
            if refs:
                for r in x:
                    if r not in refs:
                        m = "reference '{}' in resource '{}' was not recognized."
                        _LOG.warn(m.format(r, self.id))
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
    def partitioned_filepath(self):
        if self.is_abstract:
            return None
        return self.config.partitioned_dir

    @property
    def partition_source_dir(self):
        return self.normalized_filepath

    @property
    def is_abstract(self):
        return self.format is None \
               or (self.url is None and not self.url_list) \
               or self.schema is None

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

    def has_been_partitioned(self):
        part_dir_list = find_partition_dirs_for_taxonomy(self.partitioned_filepath, self.id)
        return bool(part_dir_list)

    def remove_normalize_artifacts(self):
        self._remove_taxonomy_dir(self.normalized_filepath)

    def remove_partition_artifacts(self):
        part_dir_list = find_partition_dirs_for_taxonomy(self.partitioned_filepath, self.id)
        for d in part_dir_list:
            self._remove_taxonomy_dir(d)

    def _remove_taxonomy_dir(self, directory):
        if not os.path.isdir(directory):
            return
        f_to_remove = [self.taxon_filename, 'roots.txt', 'about.json', 'details.json']
        if self.synonyms_filename:
            f_to_remove.append(self.synonyms_filename)
        _LOG.info('Removing contents of "{}"'.format(directory))
        for f in f_to_remove:
            fp = os.path.join(directory, f)
            if os.path.exists(fp):
                os.unlink(fp)
        try:
            os.rmdir(directory)
        except:
            _LOG.warn('Could not remove "{}" that directory may not be empty'.format(directory))
        else:
            _LOG.info('Removed "{}"'.format(directory))

    def get_taxdir_for_part(self, part_key):
        return get_part_inp_taxdir(self.partitioned_filepath, part_key, self.id)

    def get_taxdir_for_root_of_part(self, part_key):
        term_dir = self.get_taxdir_for_part(part_key)
        taxon_file = os.path.join(term_dir, self.taxon_filename)
        if os.path.exists(taxon_file):
            return term_dir
        par_key, misc_dir = get_par_and_par_misc_taxdir(self.partitioned_filepath,
                                                        part_key,
                                                        self.id)
        taxon_file = os.path.join(misc_dir, self.taxon_filename)
        if os.path.exists(taxon_file):
            return misc_dir
        return None

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
        schema = self.schema.lower()
        try:
            norm_fn = _schema_to_norm_fn[schema]
        except KeyError:
            m = "Normalization from \"{}\" is not currently supported"
            raise NotImplementedError(m.format(self.base_id))
        norm_fn(self.unpacked_filepath, self.normalized_filepath, self)

    def partition(self, part_name, part_keys, par_frag):
        raise NotImplementedError('Partition not implemented for base ResourceWrapper')

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
        if self.has_been_partitioned():
            part_str = "{}Has been partitioned at {}\n".format(indent, self.partitioned_filepath)
        else:
            part_str = "{}Has not been partitioned yet.\n".format(indent)
        if list_all_artifacts:
            out.write('{}{}{}{}'.format(down_str, unp_str, norm_str, part_str))
        else:
            if self.has_been_partitioned():
                out.write(part_str)
            elif self.has_been_normalized():
                out.write(norm_str)
            elif self.has_been_unpacked():
                out.write(unp_str)
            else:
                out.write(down_str)

    def get_taxon_filepath_for_part(self, fragment):
        return os.path.join(self.get_taxon_dir_for_part(fragment), self.taxon_filename)

    def get_misc_taxon_filepath_for_part(self, fragment):
        return os.path.join(self.get_misc_taxon_dir_for_part(fragment), self.taxon_filename)

    def get_taxon_dir_for_part(self, fragment):
        return get_inp_taxdir(self.partitioned_filepath, fragment, self.id)

    def get_misc_taxon_dir_for_part(self, fragment):
        return get_misc_inp_taxdir(self.partitioned_filepath, fragment, self.id)

    def get_primary_partition_map(self):
        return get_auto_gen_part_mapper(self)


# noinspection PyAbstractClass
class ExternalTaxonomyWrapper(ResourceWrapper):
    resource_type = 'external taxonomy'

    def __init__(self, obj, parent=None, refs=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs)
        # print("ET obj = {}".format(obj))

    def partition(self, part_name, part_keys, par_frag):
        do_partition(self,
                     part_name,
                     part_keys,
                     par_frag)
