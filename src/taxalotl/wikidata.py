#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .ott_schema import tax_wikidata_parser, TAXWIKIDATA_HEADER
from .taxon import Taxon
import logging

_LOG = logging.getLogger(__name__)


def _parse_taxonomy_file(taxonomy_fp):
    lp = tax_wikidata_parser
    id_2_taxon = {}
    with open(taxonomy_fp, "r") as inp:
        lit = iter(inp)
        fl = next(lit)
        assert fl == TAXWIKIDATA_HEADER
        for line in lit:
            obj = Taxon(line, line_parser=lp)
            if obj.id in id_2_taxon:
                _LOG.warning(f"Duplicate taxon ID: {obj.id}")
                if obj.__dict__ != id_2_taxon[obj.id].__dict__:
                    m = f"Duplicate taxon ID: {obj.id} with differing content."
                    raise RuntimeError(m)
            id_2_taxon[obj.id] = obj
    return id_2_taxon


def parse_wikidata(taxonomy_fp, additional_props_fp):
    id_2_taxon = _parse_taxonomy_file(taxonomy_fp)
    raise NotImplementedError("parse_wikidata")
