#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
import io
import os
import re

import csv
import sys
import logging
from enum import IntEnum

from peyutil import assure_dir_exists

from ..ott_schema import InterimTaxonomyData
from ..resource_wrapper import TaxonomyWrapper
from ..util import OutFile

_LOG = logging.getLogger(__name__)

_EXPECTED_HEADER = [
    "col:ID",
    "col:alternativeID",
    "col:nameAlternativeID",
    "col:sourceID",
    "col:parentID",
    "col:basionymID",
    "col:status",
    "col:scientificName",
    "col:authorship",
    "col:rank",
    "col:notho",
    "col:originalSpelling",
    "col:uninomial",
    "col:genericName",
    "col:infragenericEpithet",
    "col:specificEpithet",
    "col:infraspecificEpithet",
    "col:cultivarEpithet",
    "col:combinationAuthorship",
    "col:combinationAuthorshipID",
    "col:combinationExAuthorship",
    "col:combinationExAuthorshipID",
    "col:combinationAuthorshipYear",
    "col:basionymAuthorship",
    "col:basionymAuthorshipID",
    "col:basionymExAuthorship",
    "col:basionymExAuthorshipID",
    "col:basionymAuthorshipYear",
    "col:namePhrase",
    "col:nameReferenceID",
    "col:namePublishedInYear",
    "col:namePublishedInPage",
    "col:namePublishedInPageLink",
    "col:gender",
    "col:genderAgreement",
    "col:etymology",
    "col:code",
    "col:nameStatus",
    "col:accordingToID",
    "col:accordingToPage",
    "col:accordingToPageLink",
    "col:referenceID",
    "col:scrutinizer",
    "col:scrutinizerID",
    "col:scrutinizerDate",
    "col:extinct",
    "col:temporalRangeStart",
    "col:temporalRangeEnd",
    "col:environment",
    "col:species",
    "col:section",
    "col:subgenus",
    "col:genus",
    "col:subtribe",
    "col:tribe",
    "col:subfamily",
    "col:family",
    "col:superfamily",
    "col:suborder",
    "col:order",
    "col:subclass",
    "col:class",
    "col:subphylum",
    "col:phylum",
    "col:kingdom",
    "col:ordinal",
    "col:branchLength",
    "col:link",
    "col:nameRemarks",
    "col:remarks",
    "col:modified",
    "col:modifiedBy",
    "clb:merged",
]


class CoLDPIdx(IntEnum):
    ID = 0
    ALTERNATIVEID = 1
    NAMEALTERNATIVEID = 2
    SOURCEID = 3
    PARENTID = 4
    BASIONYMID = 5
    STATUS = 6
    SCIENTIFICNAME = 7
    AUTHORSHIP = 8
    RANK = 9
    NOTHO = 10
    ORIGINALSPELLING = 11
    UNINOMIAL = 12
    GENERICNAME = 13
    INFRAGENERICEPITHET = 14
    SPECIFICEPITHET = 15
    INFRASPECIFICEPITHET = 16
    CULTIVAREPITHET = 17
    COMBINATIONAUTHORSHIP = 18
    COMBINATIONAUTHORSHIPID = 19
    COMBINATIONEXAUTHORSHIP = 20
    COMBINATIONEXAUTHORSHIPID = 21
    COMBINATIONAUTHORSHIPYEAR = 22
    BASIONYMAUTHORSHIP = 23
    BASIONYMAUTHORSHIPID = 24
    BASIONYMEXAUTHORSHIP = 25
    BASIONYMEXAUTHORSHIPID = 26
    BASIONYMAUTHORSHIPYEAR = 27
    NAMEPHRASE = 28
    NAMEREFERENCEID = 29
    NAMEPUBLISHEDINYEAR = 30
    NAMEPUBLISHEDINPAGE = 31
    NAMEPUBLISHEDINPAGELINK = 32
    GENDER = 33
    GENDERAGREEMENT = 34
    ETYMOLOGY = 35
    CODE = 36
    NAMESTATUS = 37
    ACCORDINGTOID = 38
    ACCORDINGTOPAGE = 39
    ACCORDINGTOPAGELINK = 40
    REFERENCEID = 41
    SCRUTINIZER = 42
    SCRUTINIZERID = 43
    SCRUTINIZERDATE = 44
    EXTINCT = 45
    TEMPORALRANGESTART = 46
    TEMPORALRANGEEND = 47
    ENVIRONMENT = 48
    SPECIES = 49
    SECTION = 50
    SUBGENUS = 51
    GENUS = 52
    SUBTRIBE = 53
    TRIBE = 54
    SUBFAMILY = 55
    FAMILY = 56
    SUPERFAMILY = 57
    SUBORDER = 58
    ORDER = 59
    SUBCLASS = 60
    CLASS = 61
    SUBPHYLUM = 62
    PHYLUM = 63
    KINGDOM = 64
    ORDINAL = 65
    BRANCHLENGTH = 66
    LINK = 67
    NAMEREMARKS = 68
    REMARKS = 69
    MODIFIED = 70
    MODIFIEDBY = 71
    MERGED = 72


_VALID_STATUS = frozenset(
    [
        "accepted",
        "ambiguous synonym",
        "misapplied",
        "provisionally accepted",
        "synonym",
    ]
)

_ACC_STATUS = frozenset(
    [
        "accepted",
        "provisionally accepted",
    ]
)


class CDPTaxonomy(object):
    def __init__(self):
        self.by_id = {}
        self.parentless = set()
        self.to_children = {}
        self.synonyms = {}
        self.misapplied = set()

    def handle_amb_synonym(self, row):
        return self.handle_synonym(row)

    def handle_synonym(self, row):
        col_id = row[CoLDPIdx.ID]
        accepted_id = row[CoLDPIdx.PARENTID]
        name = row[CoLDPIdx.SCIENTIFICNAME]
        sset = self.synonyms.setdefault(name, set())
        sset.add((accepted_id, col_id))

    def handle_accepted(self, row):
        col_id = row[CoLDPIdx.ID]
        rank = row[CoLDPIdx.RANK]
        par_id = row[CoLDPIdx.PARENTID]
        name = row[CoLDPIdx.SCIENTIFICNAME]
        if par_id:
            self.to_children.setdefault(par_id, []).append(col_id)
        else:
            # _db_row(row)
            self.parentless.add(col_id)
        assert col_id not in self.by_id
        self.by_id[col_id] = (par_id, name, rank)

    def handle_misapplied(self, row):
        self.handle_synonym(row)
        col_id = row[CoLDPIdx.ID]
        name = row[CoLDPIdx.SCIENTIFICNAME]
        self.misapplied.add((name, col_id))

    def _write_taxon(self, outp, t_id):
        flag = ""
        par_id, name, rank = self.by_id[t_id]
        row_data = [t_id, par_id, name, rank, flag]
        row = "\t|\t".join(row_data)
        outp.write(f"{row}\n")

    def _write_taxon_tree(self, outp, root_id):
        inc = set()
        inc.add(root_id)
        self._write_taxon(outp, root_id)
        children = self.to_children.get(root_id)
        if not children:
            return inc
        for child_id in children:
            d = self._write_taxon_tree(outp, child_id)
            inc.update(d)
        return inc

    def _write_synonyms_to(self, outp, t_ids):
        for name, id_set in self.synonyms.items():
            for id_par in id_set:
                target_id, syn_id = id_par
                if target_id not in t_ids:
                    continue
                row_data = [target_id, name, ""]
                row = "\t|\t".join(row_data)
                outp.write(f"{row}\n")

    def write_to_dir(self, destination):
        inc_root = "S"
        assert inc_root in self.parentless
        vir_root = "92e52ff4-2dc6-4b35-9339-2e92035b8daf"
        assert vir_root in self.parentless

        avoid_roots = set([inc_root, vir_root])

        to_do = [i for i in self.parentless if i not in avoid_roots]
        all_roots = [to_do, [vir_root], [inc_root]]

        tax_fn = [
            "taxonomy.tsv",
            "vir_taxonomy.tsv",
            "parentless_taxonomy.tsv",
        ]
        syn_fn = [
            "synonyms.tsv",
            "vir_synonyms.tsv",
            "parentless_synonyms.tsv",
        ]
        theader = "\t|\t".join(["uid", "parent_uid", "name", "rank", "flags"])
        sheader = "\t|\t".join(["uid", "name", "type"])

        for g_idx, root_list in enumerate(all_roots):
            fn = tax_fn[g_idx]
            tax_fp = os.path.join(destination, fn)
            t_ids = set()
            # print(tax_fp)
            with open(tax_fp, "w", encoding="utf-8") as outp:
                outp.write(f"{theader}\n")
                for root in root_list:
                    n = self._write_taxon_tree(outp, root)
                    t_ids.update(n)
            fn = syn_fn[g_idx]
            syn_fp = os.path.join(destination, fn)
            with open(syn_fp, "w", encoding="utf-8") as outp:
                outp.write(f"{sheader}\n")
                self._write_synonyms_to(outp, t_ids)


def _db_row(row):
    for el in CoLDPIdx:
        print(repr(CoLDPIdx(el)), "=", repr(row[el]))


def normalize_coldp_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    taxon_fp = os.path.join(source, res_wrapper.taxon_filename)
    taxa = CDPTaxonomy()

    with open(taxon_fp, "r", encoding="utf-8") as csvfile:
        csvreader = iter(csvfile)
        hline = next(csvreader)
        assert hline[-1] == "\n"
        header = hline[:-1].split("\t")
        # If this fails, we need to make the indexing dynamic rather
        # than the CoLDPIdx enum
        assert header == _EXPECTED_HEADER

        for row_off, line in enumerate(csvreader):
            assert line[-1] == "\n"
            row = line[:-1].split("\t")
            # sys.stderr.write(f"{row_off + 1}\n")

            rstatus = row[CoLDPIdx.STATUS]
            if rstatus in _ACC_STATUS:
                taxa.handle_accepted(row)
            elif rstatus == "synonym":
                taxa.handle_synonym(row)
            elif rstatus == "ambiguous synonym":
                taxa.handle_amb_synonym(row)
            elif rstatus == "misapplied":
                taxa.handle_misapplied(row)
            else:
                assert rstatus in _VALID_STATUS
        taxa.write_to_dir(destination)
