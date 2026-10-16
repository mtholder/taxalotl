#!/usr/bin/env python
# -*- coding: utf-8 -*-
from taxalotl.resource_wrapper import TaxonomyWrapper


WIKIDATA_PARTMAP = {
    "Annelida": frozenset(["25522", "119788827"]),
    "Arachnida": frozenset(["1358"]),
    "Archaea": frozenset(["10872"]),
    "Archaeplastida": frozenset(["879246"]),
    "Arthropoda": frozenset(["1360"]),
    "Bacteria": frozenset(["10876"]),
    "Bryozoa": frozenset(["148134"]),
    "Chloroplastida": frozenset(
        [
            "132150",
        ]
    ),
    "Chordata": frozenset(["10915"]),
    "Cnidaria": frozenset(["119788837"]),
    "Coleoptera": frozenset(["22671"]),
    "Ctenophora": frozenset(["102778"]),
    "Diptera": frozenset(["25312"]),
    "Eukaryota": frozenset(["19088"]),
    "Fungi": frozenset(["764"]),
    "Glaucophyta": frozenset(["131116"]),
    "Haptophyta": frozenset(["500531"]),
    "Hymenoptera": frozenset(["22651"]),
    "Insecta": frozenset(["1390"]),
    "Life": frozenset(["10876", "10872", "Eukaryota"]),
    "Lepidoptera": frozenset(["28319"]),
    "Malacostraca": frozenset(["182978"]),
    "Metazoa": frozenset(["729"]),
    "Mollusca": frozenset(["25326"]),
    "Nematoda": frozenset(["5185"]),
    "Platyhelminthes": frozenset(["124900"]),
    "Porifera": frozenset(["18960"]),
    "Rhodophyta": frozenset(["103169"]),
    "SAR": frozenset(["137323"]),
}


class WikidataTaxonomyWrapper(TaxonomyWrapper):
    def get_primary_partition_map(self):
        if self.dynamic_part_map is not None:
            return self.dynamic_part_map
        return WIKIDATA_PARTMAP
