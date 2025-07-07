#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import tempfile
import shutil
import csv
import io
import os

from peyutil import (
    add_or_append_to_dict,
    assure_dir_exists,
    shorter_fp_form,
    write_as_json,
)
from .taxon import Taxon
from .util import OutFile
import logging

_LOG = logging.getLogger("taxalotl")

INP_OTT_TAXONOMY_HEADER = "uid\t|\tparent_uid\t|\tname\t|\trank\t|\t\n"
INP_FLAGGED_OTT_TAXONOMY_HEADER = (
    "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tflags\t|\t\n"
)
INP_FLAGGED_OTT_TAXONOMY_NO_TRAIL_HEADER = (
    "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tflags\n"
)
INP_OTT_SYNONYMS_HEADER = "uid\t|\tname\t|\ttype\t|\t\n"
FULL_OTT_HEADER = (
    "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t\n"
)
TAXWIKIDATA_HEADER = "uid\t|\tparent_uid\t|\tname\t|\trank\t|\taut_id\t|\taut_yr_id\t|\tnom_status\t|\tsrc\n"


def _parse_synonyms(tax_part):  # type (TaxonPartition) -> None
    syn_fp = tax_part.input_synonyms_filepath
    tax_part.syn_header = ""
    if not os.path.exists(syn_fp):
        return
    _LOG.debug('parsing synonyms from "{}" ...'.format(syn_fp))
    try:
        with io.open(syn_fp, "rU", encoding="utf-8") as inp:
            iinp = iter(inp)
            try:
                tax_part.syn_header = next(iinp)
            except StopIteration:
                return
            shs = tax_part.syn_header.split("\t|\t")
            if shs[0] == "uid":
                uid_ind = 0
            elif shs[1] == "uid":
                uid_ind = 1
            else:
                raise ValueError(
                    "Expected one of the first 2 columns of an OTT formatted "
                    "synonyms file to be 'uid'. Problem reading: {}".format(syn_fp)
                )
            for n, line in enumerate(iinp):
                ls = line.split("\t|\t")
                if n > 0 and n % 5000 == 0:
                    _LOG.debug(' reading synonym {:7} from "{}"'.format(n, syn_fp))
                try:
                    accept_id = ls[uid_ind]
                    try:
                        accept_id = int(accept_id)
                    except:
                        pass
                    tax_part.add_synonym(accept_id, syn_id=None, line=line)
                except:
                    _LOG.exception("Exception parsing line {}:\n{}".format(1 + n, line))
                    raise
    except:
        _LOG.exception('Exception parsing "{}"'.format(syn_fp))


def _parse_taxa(tax_part):  # type (TaxonPartition) -> None
    complete_taxon_fp = tax_part.tax_fp
    tax_part.taxon_header = ""
    if not os.path.exists(complete_taxon_fp):
        return
    ptp = shorter_fp_form(complete_taxon_fp)
    _LOG.debug('parsing taxa from "{}" ...'.format(ptp))
    with io.open(complete_taxon_fp, "r", encoding="utf-8") as inp:
        iinp = iter(inp)
        try:
            tax_part.taxon_header = next(iinp)
        except StopIteration:
            return
        for n, line in enumerate(iinp):
            ls = line.split("\t|\t")
            if n > 0 and n % 10000 == 0:
                _LOG.debug(' read taxon {:<7} from "{}" ...'.format(n, ptp))
            try:
                uid, par_id = ls[0], ls[1]
                try:
                    uid = int(uid)
                except:
                    pass
                tax_part.read_taxon_line(uid, par_id, line)
            except:
                _LOG.exception(
                    "Exception parsing line {} of {}:\n{}".format(
                        1 + n, complete_taxon_fp, line
                    )
                )
                raise


def partition_ott_by_root_id(tax_part):  # type (TaxonPartition) -> None
    _parse_synonyms(tax_part)
    _parse_taxa(tax_part)


def write_ott_taxonomy_tsv(
    out_fp,
    root_nodes,
    id_to_par,
    id_to_children,
    id_to_rank,
    id_to_name,
    id_to_flags,
    has_syn_dict,
    details_log,
    extinct_known=None,
):
    """If has_syn_dict is provided, then a list of the IDs that occur in that dict
    is returned in the order that the IDs were written to taxonomy file. This
    allows for the synonyms.tsv file to be written in a similar order, which makes browsing it
    easier.
    """
    if extinct_known is None:
        extinct_known = {}
    if has_syn_dict is None:
        has_syn_dict = {}
    syn_id_order = []
    num_tips_written = 0
    num_internals_written = 0
    rn = list(root_nodes)
    rn.sort()
    header = INP_FLAGGED_OTT_TAXONOMY_HEADER
    with OutFile(out_fp) as out:
        out.write(header)
        # need to print id, parent id, and name
        for root_id in rn:
            stack = [root_id]
            while stack:
                curr_id = stack.pop()
                if curr_id in has_syn_dict:
                    syn_id_order.append(curr_id)
                try:
                    name = id_to_name[curr_id]
                except KeyError:
                    _LOG.warning('Could not find a name for ID "{}"'.format(curr_id))
                    continue
                try:
                    par_id = id_to_par.get(curr_id)
                    if par_id is None:
                        spar_id = ""
                    else:
                        spar_id = str(par_id)
                    rank = id_to_rank.get(curr_id, "")
                    children = id_to_children.get(curr_id)
                    if children:
                        num_internals_written += 1
                        stack.extend(children)
                    else:
                        num_tips_written += 1
                    flags = id_to_flags.get(curr_id, "")
                    if flags and not isinstance(flags, str):
                        flags = ",".join(flags)
                    ev = extinct_known.get(curr_id)
                    if ev:
                        if flags:
                            flags = "{},extinct".format(flags)
                        else:
                            flags = "extinct"
                    fields = [str(curr_id), spar_id, name, rank, flags, ""]
                    try:
                        out.write("{}\n".format("\t|\t".join(fields)))
                    except:
                        _LOG.exception("error serializing {}".format(repr(fields)))
                except:
                    _LOG.error("Error writing taxon_id {}".format(curr_id))
                    raise
    details_log["num_tips_written"] = num_tips_written
    details_log["num_internals_written"] = num_internals_written
    return syn_id_order


def write_ott_synonyms_tsv(out_fp, id_to_name_name_type_list, id_order, details_log):
    num_syn_written = 0
    with OutFile(out_fp) as out:
        out.write(INP_OTT_SYNONYMS_HEADER)
        for nd_id in id_order:
            syn_list = id_to_name_name_type_list[nd_id]
            for name, name_type, syn_id in syn_list:
                num_syn_written += 1
                out.write(
                    "{}\n".format("\t|\t".join([str(nd_id), name, name_type, ""]))
                )
    details_log["num_synonyms_written"] = num_syn_written
    details_log["num_ids_with_synonyms_written"] = len(id_order)


def write_ott_forwards(out_fp, forwarded_dict):
    with OutFile(out_fp) as out:
        for key, value in forwarded_dict.items():
            out.write("{}\t{}\n".format(key, value))


def write_ncbi_details_json(fp, details_log):
    with OutFile(fp) as outs:
        write_as_json(details_log, outs, indent=2)


def read_taxonomy_to_get_id_to_name(tax_dir, id_coercion=int):
    ncbi_to_name = {}
    i = 0
    fp = os.path.join(tax_dir, "taxonomy.tsv")
    try:
        with io.open(fp, "r", encoding="utf-8") as inp:
            reader = csv.reader(inp, delimiter="\t")
            header = next(reader)
            uidx = header.index("uid")
            namex = header.index("name")
            for row in reader:
                uid = id_coercion(row[uidx])
                name = row[namex]
                if name is not None:
                    ncbi_to_name[uid] = name
                    i += 1
                    if i % 200000 == 0:
                        _LOG.info("{} {} {}".format(i, uid, name))
    except:
        _LOG.exception("error reading {}".format(fp))
        raise
    return ncbi_to_name


def int_or_str(s):
    try:
        return int(s)
    except:
        return str(s)


def raw_src_string_to_dict(raw):
    sel = raw.split(",")
    d = {}
    for el in sel:
        el = el.strip()
        if not el:
            continue
        try:
            src, sid = el.split(":")
        except ValueError:
            raise ValueError(f"Could split {el} in 2 by :")
        try:
            sid = int_or_str(sid)
        except:
            pass
        d.setdefault(src, set()).add(sid)
    return d


def full_ott_line_parser(taxon, line):
    try:
        ls = line.split("\t|\t")
        assert ls[-1] == "\n"
    except:
        _LOG.exception("Error reading line {}:\n{}".format(taxon.line_num, line))
        raise
    taxon.id = int_or_str(ls[0])
    if ls[1]:
        taxon.par_id = int_or_str(ls[1])
    else:
        taxon.par_id = None
    taxon.name = ls[2]
    if len(ls) <= 4:
        return
    if ls[3]:
        taxon.rank = ls[3]
    if len(ls) == 5:
        return
    # _LOG.debug('ls[4] = {}'.format(ls[4]))
    if ls[4]:
        d = raw_src_string_to_dict(ls[4])
        if d:
            taxon.src_dict = d
    if len(ls) == 6:
        return
    if ls[5]:
        taxon.uniqname = ls[5]
    if len(ls) > 7:
        if ls[6]:
            taxon.flags = set(ls[6].split(","))


def tax_wikidata_parser(taxon, line):
    try:
        ls = line.split("\t|\t")
        if len(ls) == 8:
            assert ls[7].endswith("\n")
            ls[7] = ls[7].strip()
        else:
            assert len(ls) > 8 and ls[-1] == "\n"
    except:
        _LOG.exception("Error reading line {}:\n{}".format(taxon.line_num, line))
        raise
    taxon.id = ls[0]
    if ls[1]:
        taxon.par_id = ls[1]
    else:
        taxon.par_id = None
    taxon.name = ls[2]
    t = taxon
    t.rank, t.aut_id, t.aut_yrs, t.nom_status, t.src_dict = None, None, None, None, None
    if ls[3]:
        taxon.rank = ls[3].strip()
    if ls[4]:
        taxon.aut_id = set(ls[4].split(","))
    if ls[5]:
        taxon.aut_yrs = set(ls[5].split(","))
    if ls[6]:
        taxon.nom_status = set(ls[6].split(","))
    d = raw_src_string_to_dict(ls[7])
    if d:
        taxon.src_dict = d


def flag_after_rank_parser(taxon, line):
    try:
        ls = line.split("\t|\t")
        if len(ls) == 5:
            assert ls[4].endswith("\n")
            ls[4] = ls[4].strip()
        else:
            assert len(ls) > 5 and ls[-1] == "\n"
    except:
        _LOG.exception("Error reading line {}:\n{}".format(taxon.line_num, line))
        raise
    taxon.id = int(ls[0])
    if ls[1]:
        taxon.par_id = int(ls[1])
    else:
        taxon.par_id = None
    taxon.name = ls[2]
    if ls[3]:
        taxon.rank = ls[3]
    if ls[4]:
        taxon.flags = set(ls[4].split(","))


HEADER_TO_LINE_PARSER = {
    FULL_OTT_HEADER: full_ott_line_parser,
    INP_OTT_TAXONOMY_HEADER: full_ott_line_parser,
    INP_FLAGGED_OTT_TAXONOMY_HEADER: flag_after_rank_parser,
    INP_FLAGGED_OTT_TAXONOMY_NO_TRAIL_HEADER: flag_after_rank_parser,
    TAXWIKIDATA_HEADER: tax_wikidata_parser,
}


# noinspection PyTypeChecker
def read_taxonomy_to_get_id_to_fields(tax_dir):
    fp = os.path.join(tax_dir, "taxonomy.tsv")
    fields = [
        "uid",
        "parent_uid",
        "name",
        "rank",
        "sourceinfo",
        "uniqname",
        "flags",
        "\n",
    ]
    expected_header = "\t|\t".join(fields)
    if not os.path.exists(fp):
        return {}
    try:
        with io.open(fp, "r", encoding="utf-8") as inp:
            iinp = iter(inp)
            header = next(iinp)
            assert header == expected_header
            id_to_obj = {}
            for n, line in enumerate(iinp):
                obj = Taxon(line, line_num=1 + n)
                oid = obj.id
                assert oid not in id_to_obj
                id_to_obj[oid] = obj
            return id_to_obj
    except:
        _LOG.exception("Error reading {}".format(fp))
        raise


def read_taxonomy_to_get_single_taxon(tax_dir, root_id):
    sri = str(root_id)
    fp = os.path.join(tax_dir, "taxonomy.tsv")
    fields = [
        "uid",
        "parent_uid",
        "name",
        "rank",
        "sourceinfo",
        "uniqname",
        "flags",
        "\n",
    ]
    expected_header = "\t|\t".join(fields)
    try:
        with io.open(fp, "r", encoding="utf-8") as inp:
            iinp = iter(inp)
            header = next(iinp)
            assert header == expected_header
            for n, line in enumerate(iinp):
                if not line.startswith(sri):
                    continue
                obj = Taxon(line, line_num=1 + n)
                if root_id == obj.id:
                    return obj
    except:
        _LOG.exception("Error reading {}".format(fp))
        raise


class InterimTaxonomyData(object):
    def __init__(self):
        self.about = {}
        self.details_log = {}
        self.forwards = {}  # from old ID to new ID
        self.to_par = {}  # ID -> parent ID or None
        self.to_children = {}  # ID to list of children IDs
        self.to_rank = {}  # ID -> rank string
        self.root_nodes = set()  # set of IDs
        self.to_name = {}  # ID -> name
        self.name_to_ids = {}  # name to
        self.to_flags = {}  # ID -> flags
        self.synonyms = {}
        self.repeated_names = set()
        self.extinct_known = None
        self.syn_id_to_valid = None
        self.extra_blob = None
        self.names_interpreted_as_changes = False

    def finalize(self):
        self.details_log["num_forwards"] = len(self.forwards)
        self.details_log["num_nodes"] = len(self.to_par)
        self.details_log["num_distinct_names"] = len(self.name_to_ids)
        self.details_log["num_ids_with_synonyms"] = len(self.synonyms)

    def register_id_and_name(self, taxon_id, name):
        self.to_name[taxon_id] = name
        if add_or_append_to_dict(self.name_to_ids, name, taxon_id):
            self.repeated_names.add(name)

    def register_synonym(self, valid_id, syn_name, name_type, syn_id=None):
        assert valid_id != syn_id
        self.synonyms.setdefault(valid_id, []).append((syn_name, name_type, syn_id))

    def fix_synonym(self, valid_id, old_valid, syn_id):
        assert valid_id != old_valid
        ol = self.synonyms.get(old_valid, [])
        matching = []
        for el in ol:
            if el[2] == syn_id:
                matching.append(el)
        nl = self.synonyms.setdefault(valid_id, [])
        for el in matching:
            ol.remove(el)
            nl.append(el)

    def write_ott_taxonomy_tsv(self, fp):
        return write_ott_taxonomy_tsv(
            fp,
            self.root_nodes,
            self.to_par,
            self.to_children,
            self.to_rank,
            self.to_name,
            self.to_flags,
            self.synonyms,
            self.details_log,
            self.extinct_known,
        )

    def write_to_dir(self, destination):
        # Write out in OTT form
        d = tempfile.mkdtemp()
        fn = [
            "taxonomy.tsv",
            "synonyms.tsv",
            "forwards.tsv",
            "about.json",
            "details.json",
        ]
        try:
            syn_order = self.write_ott_taxonomy_tsv(os.path.join(d, "taxonomy.tsv"))
            write_ott_synonyms_tsv(
                os.path.join(d, "synonyms.tsv"),
                self.synonyms,
                syn_order,
                self.details_log,
            )
            if self.forwards:
                write_ott_forwards(os.path.join(d, "forwards.tsv"), self.forwards)

            about_fp = os.path.join(d, "about.json")
            with OutFile(about_fp) as about_outs:
                write_as_json(self.about, about_outs, indent=2)
            self.finalize()
            write_ncbi_details_json(os.path.join(d, "details.json"), self.details_log)
        except:
            for f in fn:
                tf = os.path.join(d, f)
                if os.path.exists(tf):
                    try:
                        os.remove(tf)
                    except:
                        pass
            try:
                os.rmdir(d)
            except:
                pass
            raise
        assure_dir_exists(destination)
        for f in fn:
            sfp = os.path.join(d, f)
            if os.path.exists(sfp):
                dfp = os.path.join(destination, f)
                shutil.move(sfp, dfp)
        os.rmdir(d)

    def del_ids(self, id_list):
        to_name = self.to_name
        to_par = self.to_par
        to_children = self.to_children
        for taxon_id in id_list:
            if taxon_id in to_name:
                del to_name[taxon_id]
            if taxon_id in to_children:
                del to_children[taxon_id]
            pid = to_par.get(taxon_id)
            if pid:
                del to_par[taxon_id]
                pc = to_children.get(pid)
                try:
                    if pc:
                        pc.remove(taxon_id)
                except:
                    pass
