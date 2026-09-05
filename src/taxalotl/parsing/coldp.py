#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
import io
import os
import re

import csv
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
        self.synonyms.setdefault(name, set()).add((accepted_id, col_id))

    def handle_accepted(self, row):
        raise NotImplementedError("handle_accepted")

    def handle_misapplied(self, row):
        self.handle_synonym(row)
        col_id = row[CoLDPIdx.ID]
        name = row[CoLDPIdx.SCIENTIFICNAME]
        self.misapplied.add((name, col_id))


def _db_row(row):
    for el in CoLDPIdx:
        print(repr(CoLDPIdx(el)), "=", repr(row[el]))


def normalize_coldp_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    taxon_fp = os.path.join(source, res_wrapper.taxon_filename)
    taxa = CDPTaxonomy()

    with open(taxon_fp, "r", encoding="utf-8") as csvfile:
        csvreader = csv.reader(csvfile, delimiter="\t")
        header = next(csvreader)
        # If this fails, we need to make the indexing dynamic rather
        # than the CoLDPIdx enum
        assert header == _EXPECTED_HEADER

        for row in csvreader:
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

    # manifest_root = ET.parse(manifest_fp).getroot()
    # core_paths = []
    # field2index = {}
    # for el in manifest_root.findall("{http://rs.tdwg.org/dwc/text/}core"):
    #     for sub in el:
    #         if sub.tag.endswith("}id"):
    #             field2index["id"] = int(sub.attrib["index"])
    #         elif sub.tag.endswith("}field"):
    #             nns = os.path.split(sub.attrib["term"])[-1]
    #             field2index[nns] = int(sub.attrib["index"])
    #     for f in el.findall("{http://rs.tdwg.org/dwc/text/}files"):
    #         for loc in f.findall("{http://rs.tdwg.org/dwc/text/}location"):
    #             core_paths.append(loc.text.strip())
    # if len(core_paths) != 1:
    #     raise ValueError(
    #         'Did not find a single core path in DwC file ("{}") found: {}'.format(
    #             manifest_fp, core_paths
    #         )
    #     )
    # taxon_fn = core_paths[0]
    # proj_out = os.path.join(destination, "projection.tsv")
    # if not os.path.exists(proj_out):
    #     proj_in = os.path.join(source, taxon_fn)
    #     write_gbif_projection_file(proj_in, proj_out, field2index)
    # homemade = {
    #     "id": 0,
    #     "parentNameUsageID": 1,
    #     "acceptedNameUsageID": 2,
    #     "canonicalName": 3,
    #     "taxonRank": 4,
    #     "taxonomicStatus": 5,
    #     "nameAccordingTo": 6,
    # }

    # itd = InterimTaxonomyData()
    # to_remove, to_ignore, paleos = read_gbif_projection(
    #     proj_out, itd, homemade, do_gbif_checks=isinstance(res_wrapper, GBIFWrapper)
    # )
    # add_fake_root(itd)
    # remove_if_tips(itd, to_remove)
    # o_to_ignore = find_orphaned(itd)
    # to_ignore.update(o_to_ignore)
    # prune_ignored(itd, to_ignore)
    # _LOG.info("writing {} paleodb ids".format(len(paleos)))
    # with OutFile(os.path.join(destination, "paleo.tsv")) as paleofile:
    #     for taxon_id in paleos:
    #         paleofile.write("{}\n".format(taxon_id))
    # res_wrapper.post_process_interim_tax_data(itd)
    # itd.write_to_dir(destination)
