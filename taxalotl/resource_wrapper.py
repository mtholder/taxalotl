#!/usr/bin/env python
from __future__ import print_function

import io
import os
import shutil
import urllib
import urllib.request

from peyotl import (assure_dir_exists,
                    download_large_file,
                    get_logger, gunzip, gunzip_and_untar,
                    unzip)

from taxalotl.ott_schema import (INP_OTT_SYNONYMS_HEADER,
                                 INP_OTT_TAXONOMY_HEADER,
                                 partition_ott_by_root_id)
from taxalotl.newick import normalize_newick
from taxalotl.partitions import (find_partition_dirs_for_taxonomy,
                                 has_any_partition_dirs,
                                 get_auto_gen_part_mapper,
                                 get_inp_taxdir,
                                 get_misc_inp_taxdir,
                                 get_taxon_partition, )
from taxalotl.tax_partition import TAX_SLICE_CACHE

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
            content = io.open(tf, 'rU', encoding='utf-8').read()
            outfp = os.path.join(normalized_dirp, fn)
            with io.open(outfp, 'w', encoding='utf-8') as out:
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
            with io.open(tf, 'rU', encoding='utf-8') as inp:
                with io.open(outfp, 'w', encoding='utf-8') as out:
                    for line in inp:
                        ls = line.split('\t')
                        out.write('\t|\t'.join(ls))


def normalize_silva_ncbi(unpacked_dirp, normalized_dirp, resource_wrapper):
    inpfp = os.path.join(unpacked_dirp, resource_wrapper.local_filename)
    outfd = resource_wrapper.normalized_filedir
    if not os.path.exists(outfd):
        os.makedirs(outfd)
    outfp = resource_wrapper.normalized_filepath
    if resource_wrapper.schema.lower() == 'silva taxmap':
        shutil.copyfile(inpfp, outfp)
    else:
        with io.open(inpfp, 'rU', encoding='utf-8') as inp:
            with io.open(outfp, 'w', encoding='utf-8') as outp:
                for line in inp:
                    ls = line.strip()
                    if not ls:
                        continue
                    if ls[0] != '>':
                        raise ValueError(
                            'Expecting each line to start with > found:\n{}'.format(ls))
                    spl = ls.split(' ')
                    silva_info = spl[0][1:]
                    rest = ' '.join(spl[1:])
                    sispl = silva_info.split('.')
                    first_col = '.'.join(sispl[:-2])
                    sec_col, third_col = sispl[-2:]
                    label = rest.split(';')[-1]
                    if not rest.endswith(';'):
                        rest += ';'
                    ol = '\t'.join([first_col, sec_col, third_col, rest, label])
                    outp.write('{}\n'.format(ol))


_schema_to_norm_fn = {"headerless ott": copy_and_add_ott_headers,
                      "newick": normalize_newick,
                      "ott": copy_taxonomy_by_linking,
                      "ott id csv": copy_id_list_by_linking,
                      "silva_ncbi": normalize_silva_ncbi,
                      'fasta silva taxmap': normalize_silva_ncbi,
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
    _norm_filename = taxon_filename
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

    def get_self_and_children(self):
        # type: () -> [ResourceWrapper]
        r = [self]
        if self.children:
            for c in self.children:
                r.extend(c.get_self_and_children())
        return r

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
    def normalized_filedir(self):
        if self.is_abstract:
            return None
        return os.path.join(self.config.normalized_dir, self.id)

    @property
    def normalized_filepath(self):
        fd = self.normalized_filedir
        if fd is None:
            return None
        return os.path.join(fd, self._norm_filename)

    @property
    def partitioned_filepath(self):
        if self.is_abstract:
            return None
        return self.config.partitioned_dir

    @property
    def partition_source_dir(self):
        return self.normalized_filedir

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
        dfp = self.normalized_filedir
        return (dfp is not None
                and os.path.exists(dfp)
                and os.path.exists(self.normalized_filepath))

    def has_been_partitioned(self):
        return has_any_partition_dirs(self.partitioned_filepath, self.id)

    def remove_normalize_artifacts(self):
        self._remove_taxonomy_dir(self.normalized_filedir)

    def remove_partition_artifacts(self):
        part_dir_list = find_partition_dirs_for_taxonomy(self.partitioned_filepath, self.id)
        for d in part_dir_list:
            self._remove_taxonomy_dir(d)

    def _remove_taxonomy_dir(self, directory):
        if not os.path.isdir(directory):
            return
        f_to_remove = [self.taxon_filename,
                       '__roots__.json',
                       'about.json',
                       'details.json',
                       '__accum_des__.json',
                       ]
        if self.synonyms_filename:
            f_to_remove.append(self.synonyms_filename)
        for f in f_to_remove:
            fp = os.path.join(directory, f)
            if os.path.exists(fp):
                _LOG.info('Removing "{}"'.format(fp))
                os.unlink(fp)
        try:
            os.rmdir(directory)
        except:
            _LOG.warn('Could not remove "{}" that directory may not be empty'.format(directory))
        else:
            _LOG.info('Removed "{}"'.format(directory))

    def get_taxdir_for_part(self, part_key):
        return self.config.get_part_inp_taxdir(part_key, self.id)

    def get_taxdir_for_root_of_part(self, part_key):
        term_dir = self.get_taxdir_for_part(part_key)
        taxon_file = os.path.join(term_dir, self.taxon_filename)
        if os.path.exists(taxon_file):
            return term_dir
        par_key, misc_dir = self.config.get_par_and_par_misc_taxdir(part_key, self.id)
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
            with urllib.request.urlopen(self.url) as req:
                with io.open(dfp, 'w', 'utf-8') as outp:
                    outp.write(req.read())
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
        norm_fn(self.unpacked_filepath, self.normalized_filedir, self)

    def write_status(self, out, indent='', list_all_artifacts=False, hanging_indent='  '):
        dfp = self.download_filepath
        src_str = '{} '.format(self.source) if self.source else ''
        if dfp is None:
            lo = self.get_leaf_obj()
            if lo == self:
                suff = ''
            else:
                suff = ' the most recent version appears to be "{}"'.format(lo.id)
            template = "{}{}: {}. {}This is a class of resource, not a downloadable artifact{}.\n"
            out.write(template.format(indent, self.id, self.resource_type, src_str, suff))
            return
        out.write('{}{}: {} {}'.format(indent, self.id, self.resource_type, src_str))
        if self.version:
            out.write('version {}. '.format(self.version))
        else:
            out.write('(unversioned). ')
        out.write("date={}\n".format(self.date if self.date else 'unknown'))
        hi = '{}{}'.format(indent, hanging_indent)
        s = "is at" if self.has_been_downloaded() else "not yet downloaded to"
        down_str = "{}Raw ({} format) {} {}\n".format(hi, self.format, s, dfp)
        ufp = self.unpacked_filepath
        s = "is at" if self.has_been_unpacked() else "not yet unpacked to"
        unp_str = "{}Raw ({} schema) {} {}\n".format(hi, self.schema, s, ufp)
        nfp = self.normalized_filedir
        s = "is at" if self.has_been_normalized() else "not yet normalized to"
        norm_str = "{}OTT formatted form {} {}\n".format(hi, s, nfp)
        if self.has_been_partitioned():
            part_str = "{}Has been partitioned at {}\n".format(hi, self.partitioned_filepath)
        else:
            part_str = "{}Has not been partitioned yet.\n".format(hi)
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

    def has_been_partitioned_for_fragment(self, fragment):
        return os.path.exists(self.get_misc_taxon_filepath_for_part(fragment))

    @property
    def base_resource(self):
        return self if not self.parent else self.parent.base_resource

    @property
    def alias_list(self):
        x = getattr(self, 'aliases', [])
        return list(x) if x else []


class AbstractResourceWrapper(ResourceWrapper):
    def __init__(self, obj, parent=None, refs=None, config=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs, config=config)


# noinspection PyAbstractClass
class TaxonomyWrapper(ResourceWrapper):
    resource_type = 'external taxonomy'
    schema = {'headerless ott', "newick", "ott", "ott id csv", "tab-separated ott"}

    def __init__(self, obj, parent=None, refs=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs)
        # print("ET obj = {}".format(obj))

    def accumulate_separated_descendants(self, scaffold_dir):
        scaffold_anc = os.path.split(scaffold_dir)[0]
        pd = self.partitioned_filepath
        # _LOG.debug('comparing "{}" and "{}"'.format(pd, scaffold_anc))
        if scaffold_anc == pd:
            return
        pd = '{}/'.format(pd)
        assert scaffold_anc.startswith(pd)
        frag = scaffold_dir[len(pd):]
        _LOG.info('frag = {}'.format(frag))
        tax_part = get_taxon_partition(self, fragment=frag)
        tax_part.read_inputs_for_read_only()
        # root_ids = tax_part.get_root_ids()
        forest = tax_part.get_taxa_as_forest()
        des_accum = tax_part.read_acccumulated_des()
        if des_accum:
            for des_id, par_id, line in des_accum:
                assert forest.get_taxon(par_id) is not None
        accum_list = []
        pass_through = tax_part.read_pass_through_des()
        if pass_through:
            accum_list.extend(pass_through)
        crs = set()
        for tree in forest.trees:
            root = tree.root
            crs.add(root.id)
            accum_list.append((root.id, root.par_id, root.line))
        # _LOG.info('accum_list: "{}"'.format(accum_list))
        anc_frag = os.path.split(frag)[0]
        while not self.has_been_partitioned_for_fragment(anc_frag):
            if not anc_frag:
                assert False
            anc_frag = os.path.split(anc_frag)[0]
        anc_tax_part = get_taxon_partition(self, anc_frag)
        anc_tax_part.register_accumulated_des(accum_list)
        TAX_SLICE_CACHE.try_del(tax_part.cache_key)

    @property
    def is_abstract_input_resource_type(self):
        return self.id == self.base_id and self.base_id != 'ott'

    @property
    def unversioned_base_name(self):
        if self.base_id == 'irmng_ot':
            return 'irmng'  # sorry this is hacky...
        return self.base_id