#!/usr/bin/env python
from __future__ import print_function

import io
import os
import shutil
import urllib
import urllib.request
import logging

from peyutil import (
    assure_dir_exists,
    download_large_file,
    gunzip,
    gunzip_and_untar,
    unzip,
    read_as_json,
)

from .ott_schema import (
    INP_OTT_SYNONYMS_HEADER,
    INP_OTT_TAXONOMY_HEADER,
    partition_ott_by_root_id,
    INP_FLAGGED_OTT_TAXONOMY_HEADER,
)
from .newick import normalize_newick
from .cmds.partitions import (
    find_partition_dirs_for_taxonomy,
    get_all_partition_dirs,
    has_any_partition_dirs,
    get_inp_taxdir,
    get_misc_inp_taxdir,
)
from .tax_partition import TAX_SLICE_CACHE, ROOTS_FILENAME, ACCUM_DES_FILENAME
from .util import unlink, OutFile, OutDir
from .wikispecies import parse_wikispecies
from .wikidata import parse_wikidata

_LOG = logging.getLogger(__name__)


######################################################################################
def unpack_archive(archive_fp, unpack_fp, archive_format, wrapper):
    afl = archive_format.lower()
    if afl in ["tar+gzip"]:
        _LOG.debug("gunzip_and_untar from {} to {} ...".format(archive_fp, unpack_fp))
        gunzip_and_untar(archive_fp, unpack_fp)
        _LOG.debug("gunzip_and_untar from {} to {} done.".format(archive_fp, unpack_fp))
    elif afl == "zip":
        _LOG.debug("unzip from {} to {} ...".format(archive_fp, unpack_fp))
        unzip(archive_fp, unpack_fp)
        _LOG.debug("unzip from {} to {} done.".format(archive_fp, unpack_fp))
    elif afl == "gzip":
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
    elif afl == "text":
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
OTT_TAXONOMY_FILENAMES = (
    "about.json",
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
    "version.txt",
)

OTT_TAXONOMY_ID_FILES = ("by_qid.csv",)


def copy_file_list_by_linking(unpacked_dirp, normalized_dirp, file_list):
    assure_dir_exists(normalized_dirp)
    for fn in file_list:
        ufp = os.path.join(unpacked_dirp, fn)
        if os.path.exists(ufp):
            dfp = os.path.join(normalized_dirp, fn)
            if os.path.exists(dfp):
                _LOG.info(
                    'File already exists at "{}". Skipping link creation.'.format(dfp)
                )
            else:
                os.symlink(ufp, dfp)


# noinspection PyUnusedLocal
def copy_taxonomy_by_linking(unpacked_dirp, normalized_dirp, resource_wrapper):
    copy_file_list_by_linking(unpacked_dirp, normalized_dirp, OTT_TAXONOMY_FILENAMES)


# noinspection PyUnusedLocal
def copy_id_list_by_linking(unpacked_dirp, normalized_dirp, resource_wrapper):
    copy_file_list_by_linking(unpacked_dirp, normalized_dirp, OTT_TAXONOMY_ID_FILES)


# noinspection PyUnusedLocal
def copy_and_add_ott_headers(unpacked_dirp, normalized_dirp, resource_wrapper):
    motf = list(OTT_TAXONOMY_FILENAMES)
    special = [
        ("taxonomy.tsv", INP_OTT_TAXONOMY_HEADER),
        ("synonyms.tsv", INP_OTT_SYNONYMS_HEADER),
    ]
    for fn, header in special:
        motf.remove(fn)
    copy_file_list_by_linking(unpacked_dirp, normalized_dirp, motf)
    for fn, header in special:
        tf = os.path.join(unpacked_dirp, fn)
        if os.path.isfile(tf):
            content = io.open(tf, "r", encoding="utf-8").read()
            outfp = os.path.join(normalized_dirp, fn)
            with OutFile(outfp) as out:
                out.write(header)
                out.write(content)


# noinspection PyUnusedLocal
def normalize_tab_sep_ott(unpacked_dirp, normalized_dirp, resource_wrapper):
    motf = list(OTT_TAXONOMY_FILENAMES)
    special = [("taxonomy.tsv", INP_OTT_TAXONOMY_HEADER)]
    for fn, header in special:
        motf.remove(fn)
    copy_file_list_by_linking(unpacked_dirp, normalized_dirp, motf)
    for fn, header in special:
        tf = os.path.join(unpacked_dirp, fn)
        if os.path.isfile(tf):
            outfp = os.path.join(normalized_dirp, fn)
            with io.open(tf, "r", encoding="utf-8") as inp:
                with OutFile(outfp) as out:
                    for line in inp:
                        ls = line.split("\t")
                        out.write("\t|\t".join(ls))


def normalize_wikispecies(unpacked_dirp, normalized_dirp, resource_wrapper):
    # import sys; sys.exit(f"{unpacked_dirp}\n{normalized_dirp}\n{resource_wrapper.__dict__}\n")
    fn = "pages-meta-current.xml"
    infp = os.path.join(unpacked_dirp, fn)
    if not os.path.isfile(infp):
        raise RuntimeError(f"{infp} does not exist")
    parse_wikispecies(infp)


def normalize_wikidata(unpacked_dirp, normalized_dirp, resource_wrapper):
    # import sys; sys.exit(f"{unpacked_dirp}\n{normalized_dirp}\n{resource_wrapper.__dict__}\n")
    fn = "taxonomy.tsv"
    infp = os.path.join(unpacked_dirp, fn)
    if not os.path.isfile(infp):
        raise RuntimeError(f"{infp} does not exist")
    add_fp = os.path.join(unpacked_dirp, "additional-properties.tsv")
    id_2_taxon, synonyms = parse_wikidata(infp, add_fp)
    for syn_id in synonyms:
        syn_taxon = id_2_taxon[syn_id]
        valid_syn = set()
        for pvi in syn_taxon.synonyms:
            if pvi not in synonyms:
                valid_syn.add(pvi)
        if len(valid_syn) == 0:
            m = f"Taxon {syn_id} has synonym(s): {syn_taxon.synonyms}, but size of valid set ({valid_syn}) is not 1"
            _LOG.warning(m)
        else:
            syn_taxon.valid_syn = valid_syn
    outfd = resource_wrapper.normalized_filedir
    with OutDir(outfd):
        pass
    outfp = resource_wrapper.normalized_filepath
    with OutFile(outfp) as out:
        out.write(INP_FLAGGED_OTT_TAXONOMY_HEADER)
        for taxon in id_2_taxon.values():
            if hasattr(taxon, "valid_syn"):
                continue
            pid = taxon.par_id[1:].strip() if taxon.par_id else ""
            r = taxon.rank if taxon.rank else ""
            els = [taxon.id[1:].strip(), pid, taxon.name, r, taxon.get_flag_str()]
            row = "\t|\t".join(els)
            out.write(f"{row}\n")
    syn_fp = os.path.join(outfd, resource_wrapper.synonyms_filename)
    with OutFile(syn_fp) as out:
        out.write(INP_OTT_SYNONYMS_HEADER)
        for syn in synonyms:
            taxon = id_2_taxon[syn]
            for referred_id in taxon.synonyms:
                try:
                    ref_tax = id_2_taxon[referred_id]
                except KeyError:
                    continue
                if hasattr(ref_tax, "valid_syn"):
                    continue
                num_id = referred_id[1:].strip()
                row = "\t|\t".join([num_id, taxon.name, ""])
                out.write(f"{row}\n")


def normalize_silva_ncbi(unpacked_dirp, normalized_dirp, resource_wrapper):
    inpfp = os.path.join(unpacked_dirp, resource_wrapper.local_filename)
    outfd = resource_wrapper.normalized_filedir
    with OutDir(outfd):
        pass
    outfp = resource_wrapper.normalized_filepath
    if resource_wrapper.schema.lower() == "silva taxmap":
        shutil.copyfile(inpfp, outfp)
    else:
        with io.open(inpfp, "r", encoding="utf-8") as inp:
            with OutFile(outfp) as outp:
                for line in inp:
                    ls = line.strip()
                    if not ls:
                        continue
                    if ls[0] != ">":
                        raise ValueError(
                            "Expecting each line to start with > found:\n{}".format(ls)
                        )
                    spl = ls.split(" ")
                    silva_info = spl[0][1:]
                    rest = " ".join(spl[1:])
                    sispl = silva_info.split(".")
                    first_col = ".".join(sispl[:-2])
                    sec_col, third_col = sispl[-2:]
                    label = rest.split(";")[-1]
                    if not rest.endswith(";"):
                        rest += ";"
                    ol = "\t".join([first_col, sec_col, third_col, rest, label])
                    outp.write("{}\n".format(ol))


_schema_to_norm_fn = {
    "headerless ott": copy_and_add_ott_headers,
    "newick": normalize_newick,
    "ott": copy_taxonomy_by_linking,
    "ott id csv": copy_id_list_by_linking,
    "silva_ncbi": normalize_silva_ncbi,
    "fasta silva taxmap": normalize_silva_ncbi,
    "tab-separated ott": normalize_tab_sep_ott,
    "wikispecies": normalize_wikispecies,
    "taxalotlwikidata": normalize_wikidata,
}

_known_res_attr = frozenset(
    [
        "aliases",
        "base_id",  # base ID in inherits from graph
        "copy_status",
        "date",
        "depends_on",
        "doc_url",
        "format",
        "id",
        "inherits_from",
        "inputs",
        "id_list",
        "latest_download_url",
        "license_url",
        "license_or_tou_info",
        "local_filename",
        "maintainer",
        "notes",
        "preceded_by",
        "references",
        "resource_type",
        "schema",
        "source",
        "stats",
        "url",
        "url_list",
        "version",
    ]
)


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
    taxon_filename = "taxonomy.tsv"
    _norm_filename = taxon_filename
    synonyms_filename = "synonyms.tsv"
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
        self.id = obj["id"]
        if self.references:
            x = (
                [self.references]
                if not isinstance(self.references, list)
                else self.references
            )
            if refs:
                for r in x:
                    if r not in refs:
                        m = "reference '{}' in resource '{}' was not recognized."
                        _LOG.warning(m.format(r, self.id))
        self._config = config

    @property
    def config(self):
        if self._config is None:
            m = (
                "Inadequately initialized ResourceWrapper: "
                "config attribute of a resource must be set."
            )
            raise RuntimeError(m)
        return self._config

    @config.setter
    def config(self, c):
        self._config = c

    def get_part_filepaths_by_name(self, basename):
        """Returns a list filepaths ending in /basename for this partitioned res"""
        part_dirs = self.get_partitions_roots()
        rfp_list = []
        for opd in part_dirs:
            rfp = os.path.join(opd, basename)
            if os.path.isfile(rfp):
                rfp_list.append(rfp)
        return rfp_list

    def get_part_roots_filepaths(self):
        """Returns a list of all of the __roots__.json files partitioned."""
        return self.get_part_filepaths_by_name(ROOTS_FILENAME)

    def get_part_taxa_filepaths(self):
        """Returns a list of all of the "taxonomy.tsv" files partitioned."""
        return self.get_part_filepaths_by_name("taxonomy.tsv")

    def get_part_syn_filepaths(self):
        """Returns a list of all of the "synonyms.tsv" files partitioned."""
        return self.get_part_filepaths_by_name("synonyms.tsv")

    def get_part_clade_names_and_blobs(self):
        by_name = {}
        rfp_list = self.get_part_roots_filepaths()
        for rfp in rfp_list:
            rfp = os.path.join(opd, ROOTS_FILENAME)
            blob = read_as_json(rfp)
            if len(blob) != 1:
                raise RuntimeError(f"Multiple roots found at {blob}")
            val = [i for i in blob.values()][0]
            vn = val["name"]
            by_name[vn] = val
        return by_name

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
        return (
            self.format is None
            or (self.url is None and not self.url_list)
            or self.schema is None
        )

    @property
    def is_abstract_input_resource_type(self):
        return False

    def has_been_downloaded(self):
        dfp = self.download_filepath
        return dfp is not None and os.path.exists(dfp)

    def has_been_unpacked(self):
        dfp = self.unpacked_filepath
        return dfp is not None and os.path.exists(dfp)

    def has_been_normalized(self):
        dfp = self.normalized_filedir
        return (
            dfp is not None
            and os.path.exists(dfp)
            and os.path.exists(self.normalized_filepath)
        )

    def has_been_partitioned(self):
        return has_any_partition_dirs(self.partitioned_filepath, self.id)

    def get_partitions_roots(self):
        return get_all_partition_dirs(self.partitioned_filepath, self.id)

    def remove_normalize_artifacts(self):
        self._remove_taxonomy_dir(self.normalized_filedir)

    def remove_partition_artifacts(self):
        part_dir_list = find_partition_dirs_for_taxonomy(
            self.partitioned_filepath, self.id
        )
        for d in part_dir_list:
            self._remove_taxonomy_dir(d)

    def _remove_taxonomy_dir(self, directory):
        if not os.path.isdir(directory):
            return
        f_to_remove = [
            self.taxon_filename,
            "taxonomy.tsv",
            ROOTS_FILENAME,
            "about.json",
            "details.json",
            ACCUM_DES_FILENAME,
        ]
        if self.synonyms_filename:
            f_to_remove.append(self.synonyms_filename)
        for f in f_to_remove:
            fp = os.path.join(directory, f)
            if os.path.exists(fp):
                unlink(fp)
        try:
            os.rmdir(directory)
        except:
            _LOG.warning(
                f'Could not remove "{directory}" that directory may not be empty'
            )
        else:
            _LOG.info('Removed "{}"'.format(directory))

    def download(self):
        dfp = self.download_filepath
        if dfp is None:
            m = "Resource {} appears to be abstract, therefore not downloadable"
            raise RuntimeError(m.format(self.id))
        if self.url.startswith("ftp://"):
            _LOG.debug("Starting FTP download from {} to {}".format(self.url, dfp))
            with urllib.request.urlopen(self.url) as req:
                with OutFile(dfp, mode="wb") as outp:
                    outp.write(req.read())
            _LOG.debug("Download from {} to {} completed.".format(self.url, dfp))
        else:
            _LOG.debug("Starting download from {} to {}".format(self.url, dfp))
            download_large_file(self.url, dfp)
            _LOG.debug("Download from {} to {} completed.".format(self.url, dfp))

    def unpack(self):
        unpack_archive(
            self.download_filepath, self.unpacked_filepath, self.format, self
        )

    def normalize(self):
        schema = self.schema.lower()
        try:
            norm_fn = _schema_to_norm_fn[schema]
        except KeyError:
            m = f'Normalization from "{self.base_id}" ({schema}) is not currently supported'
            raise NotImplementedError(m)
        norm_fn(self.unpacked_filepath, self.normalized_filedir, self)

    def write_status(
        self, out, indent="", list_all_artifacts=False, hanging_indent="  "
    ):
        dfp = self.download_filepath
        src_str = "{} ".format(self.source) if self.source else ""
        if dfp is None:
            lo = self.get_leaf_obj()
            if lo == self:
                suff = ""
            else:
                suff = ' the most recent version appears to be "{}"'.format(lo.id)
            template = "{}{}: {}. {}This is a class of resource, not a downloadable artifact{}.\n"
            out.write(
                template.format(indent, self.id, self.resource_type, src_str, suff)
            )
            return
        out.write("{}{}: {} {}".format(indent, self.id, self.resource_type, src_str))
        if self.version:
            out.write("version {}. ".format(self.version))
        else:
            out.write("(unversioned). ")
        out.write("date={}\n".format(self.date if self.date else "unknown"))
        hi = "{}{}".format(indent, hanging_indent)
        s = "is at" if self.has_been_downloaded() else "not yet downloaded to"
        down_str = "{}Raw ({} format) {} {}\n".format(hi, self.format, s, dfp)
        ufp = self.unpacked_filepath
        s = "is at" if self.has_been_unpacked() else "not yet unpacked to"
        unp_str = "{}Raw ({} schema) {} {}\n".format(hi, self.schema, s, ufp)
        nfp = self.normalized_filedir
        s = "is at" if self.has_been_normalized() else "not yet normalized to"
        norm_str = "{}OTT formatted form {} {}\n".format(hi, s, nfp)
        if self.has_been_partitioned():
            part_str = "{}Has been partitioned at {}\n".format(
                hi, self.partitioned_filepath
            )
        else:
            part_str = "{}Has not been partitioned yet.\n".format(hi)
        if list_all_artifacts:
            out.write("{}{}{}{}".format(down_str, unp_str, norm_str, part_str))
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
        return os.path.join(self.get_taxon_dir_for_part(fragment), "taxonomy.tsv")

    def get_misc_taxon_filepath_for_part(self, fragment):
        return os.path.join(self.get_misc_taxon_dir_for_part(fragment), "taxonomy.tsv")

    def get_taxon_dir_for_part(self, fragment):
        return get_inp_taxdir(self.partitioned_filepath, fragment, self.id)

    def get_misc_taxon_dir_for_part(self, fragment):
        return get_misc_inp_taxdir(self.partitioned_filepath, fragment, self.id)

    def get_primary_partition_map(self):
        raise RuntimeError(
            f"Resource {self.id} does not contain a hard-coded primary partition map"
        )

    def has_been_partitioned_for_fragment(self, fragment):
        return os.path.exists(self.get_misc_taxon_filepath_for_part(fragment))

    @property
    def base_resource(self):
        return self if not self.parent else self.parent.base_resource

    @property
    def alias_list(self):
        x = getattr(self, "aliases", [])
        return list(x) if x else []


class AbstractResourceWrapper(ResourceWrapper):
    def __init__(self, obj, parent=None, refs=None, config=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs, config=config)

    @property
    def is_abstract_input_resource_type(self):
        return True


# noinspection PyAbstractClass
class TaxonomyWrapper(ResourceWrapper):
    resource_type = "external taxonomy"
    schema = {"headerless ott", "newick", "ott", "ott id csv", "tab-separated ott"}

    def __init__(self, obj, parent=None, refs=None, config=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs, config=config)
        self.part_name_to_tax_part_in_mem = {}
        # print("ET obj = {}".format(obj))

    @property
    def is_abstract_input_resource_type(self):
        return self.id == self.base_id and self.base_id != "ott"

    @property
    def unversioned_base_name(self):
        if self.base_id == "irmng_ot":
            return "irmng"  # sorry this is hacky...
        return self.base_id

    def post_process_interim_tax_data(self, interim_tax_data):
        pass

    def collapse_as_incertae_sedis_interim_tax_data(self, interim_tax_data, prefix):
        if interim_tax_data.names_interpreted_as_changes:
            return
        container_id_to_par_id = {}
        for old_id, name in interim_tax_data.to_name.items():
            if name.lower().startswith(prefix):
                container_id_to_par_id[old_id] = interim_tax_data.to_par[old_id]
        if container_id_to_par_id:
            # deal with possibility of nested containters
            while True:
                new_cont_map = {}
                for k, v in container_id_to_par_id.items():
                    if v in container_id_to_par_id:
                        new_cont_map[k] = container_id_to_par_id[v]
                if new_cont_map:
                    container_id_to_par_id.update(new_cont_map)
                else:
                    break

            new_par_id = {}
            for old_id, par_id in interim_tax_data.to_par.items():
                np = container_id_to_par_id.get(par_id)
                if np:
                    new_par_id[old_id] = np
            for old_id, new_par in new_par_id.items():
                interim_tax_data.to_children[new_par].append(old_id)
                interim_tax_data.to_par[old_id] = new_par
                interim_tax_data.to_flags.setdefault(old_id, []).append(
                    "incertae_sedis"
                )
            interim_tax_data.del_ids(container_id_to_par_id.keys())
        interim_tax_data.names_interpreted_as_changes = True


class GenericTaxonomyWrapper(TaxonomyWrapper):
    def __init__(self, res_id, config):
        super(GenericTaxonomyWrapper, self).__init__(
            {"id": res_id}, None, config=config
        )

    @property
    def is_abstract(self):
        return False
