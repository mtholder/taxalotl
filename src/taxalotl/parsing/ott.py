#!/usr/bin/env python
# from __future__ import print_function

import os
import logging

from peyutil import StringIO

from ..cmds.partitions import (
    MISC_DIRNAME,
    PREORDER_PART_LIST,
    NAME_TO_PARTS_SUBSETS,
)
from ..resource_wrapper import ResourceWrapper, TaxonomyWrapper

_LOG = logging.getLogger(__name__)

OTT_PARTMAP = {
    "Archaea": frozenset(["996421"]),
    "Bacteria": frozenset(["844192"]),
    "Eukaryota": frozenset(["304358"]),
    "SAR": frozenset(["5246039"]),
    "Haptophyta": frozenset(["151014"]),
    "Rhodophyta": frozenset(["878953"]),
    "Archaeplastida": frozenset(["5268475"]),
    "Glaucophyta": frozenset(["664970"]),
    "Chloroplastida": frozenset(["361838"]),
    "Fungi": frozenset(["352914"]),
    "Metazoa": frozenset(["691846"]),
    "Annelida": frozenset(["941620"]),
    "Arthropoda": frozenset(["632179"]),
    "Malacostraca": frozenset(["212701"]),
    "Arachnida": frozenset(["511967"]),
    "Insecta": frozenset(["1062253"]),
    "Diptera": frozenset(["661378"]),
    "Coleoptera": frozenset(["865243"]),
    "Lepidoptera": frozenset(["965954"]),
    "Hymenoptera": frozenset(["753726"]),
    "Bryozoa": frozenset(["442934"]),
    "Chordata": frozenset(["125642"]),
    "Cnidaria": frozenset(["641033"]),
    "Ctenophora": frozenset(["641212"]),
    "Mollusca": frozenset(["802117"]),
    "Nematoda": frozenset(["395057"]),
    "Platyhelminthes": frozenset(["555379"]),
    "Porifera": frozenset(["67819"]),
    "Viruses": frozenset(["4807313"]),
}

# Unused separation taxa: cellular organisms	93302


OTT_3_SEPARATION_TAXA = OTT_PARTMAP


UNSTABLE_SRC_PREFIXES = frozenset(["h2007", "study713", "https", "http"])


class NewSeparator(object):
    def __init__(self, ott_taxon_obj):
        self.taxon = ott_taxon_obj
        self.sub_separators = {}

    def __str__(self):
        out = StringIO()
        self.write_str(out)
        return out.getvalue()

    def write_str(self, out, pref=""):
        t = self.taxon
        if pref:
            pref = pref + "/"
        out.write(
            "{}{} ott{} sources={}\n".format(
                pref, t.name_that_is_unique, t.id, t.src_dict
            )
        )
        for name, sep in self.sub_separators.items():
            if pref:
                ns = pref + name
            else:
                ns = name
            sep.write_str(out, ns)

    def as_dict(self):
        d = {"name": self.taxon.name, "uniqname": self.taxon.name_that_is_unique}
        sd = {}
        for k, v in self.taxon.src_dict.items():
            vl = list(v)
            vl.sort()
            sd[k] = vl
        d["src_dict"] = sd
        s = {}
        for v in self.sub_separators.values():
            s[v.taxon.id] = v.as_dict()
        if s:
            d["sub"] = s
        return d

    def num_sub_separators(self):
        n = 0
        for el in self.sub_separators.values():
            n += 1 + el.num_sub_separators()
        return n


def _add_nst_subtree_el_to_dict(rd, nst_el, par_to_child):
    sep_taxon, children = nst_el
    if sep_taxon is not None:
        nst = NewSeparator(sep_taxon)
        nd = nst.sub_separators
        rd[sep_taxon.name_that_is_unique] = nst
    else:
        nd = rd
    for c in children:
        next_el = par_to_child[c]
        _add_nst_subtree_el_to_dict(nd, next_el, par_to_child)


NON_SEP_RANKS = frozenset(
    [
        "forma",
        "no rank - terminal",
        "species",
        "species group",
        "species subgroup",
        "varietas",
        "variety",
    ]
)
MIN_SEP_SIZE = 20


DEFAULT_REL_SRC_SET = frozenset(["gbif", "irmng", "ncbi", "worms"])
NO_WORMS_REL_SRC_SET = frozenset(["gbif", "irmng", "ncbi"])
SILVA_NOT_WORMS_REL_SRC_SET = frozenset(["gbif", "irmng", "ncbi", "silva"])
PART_KEY_TO_REL_SRC_SET = {
    "Insecta": NO_WORMS_REL_SRC_SET,
    "Hymenoptera": NO_WORMS_REL_SRC_SET,
    "Diptera": NO_WORMS_REL_SRC_SET,
    "Coleoptera": NO_WORMS_REL_SRC_SET,
    "Lepidoptera": NO_WORMS_REL_SRC_SET,
    "Chordata": NO_WORMS_REL_SRC_SET,
    "Fungi": NO_WORMS_REL_SRC_SET,
    "Bacteria": SILVA_NOT_WORMS_REL_SRC_SET,
    "Archaea": SILVA_NOT_WORMS_REL_SRC_SET,
}


# noinspection PyAbstractClass
class OTTaxonomyWrapper(TaxonomyWrapper):
    resource_type = "open tree taxonomy"
    schema = {"ott"}

    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    def get_primary_partition_map(self):
        return OTT_PARTMAP


# noinspection PyAbstractClass
class OTTaxonomyIdListWrapper(ResourceWrapper):
    resource_type = "open tree taxonomy idlist"
    schema = {resource_type, "ott id csv"}
    _norm_filename = "by_qid.csv"

    def __init__(self, obj, parent=None, refs=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs)
