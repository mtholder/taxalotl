#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .ott_schema import tax_wikidata_parser, TAXWIKIDATA_HEADER
from .taxon import Taxon


def _parse_taxonomy_file(taxonomy_fp):
    lp = tax_wikidata_parser
    taxa = []
    with open(taxonomy_fp, "r") as inp:
        lit = iter(inp)
        fl = next(lit)
        assert fl == TAXWIKIDATA_HEADER
        for line in lit:
            obj = Taxon(line, line_parser=lp)
            taxa.append(obj)
    return taxa


def parse_wikidata(taxonomy_fp, additional_props_fp):
    taxa = _parse_taxonomy_file(taxonomy_fp)
    raise NotImplementedError("parse_wikidata")
