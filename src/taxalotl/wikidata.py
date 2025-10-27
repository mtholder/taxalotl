#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .ott_schema import tax_wikidata_parser, TAXWIKIDATA_HEADER
from .taxon import Taxon
import logging
from qwikidata.entity import WikidataItem


_LOG = logging.getLogger(__name__)

P_IS_INSTANCE_OF = "P31"
Q_TAXON = "Q16521"
Q_FOSSIL_TAXON = "Q23038290"
Q_MONOTYPIC_FOSSIL_TAXON = "Q47487597"
Q_CLADE = "Q713623"
Q_MONOTYPIC_TAXON = "Q310890"
Q_SPECIES = "Q15458879"
Q_SUBSPECIES = "Q68947"
Q_MONOPHYLY = "Q210958"  # Sigh...
Q_LAZARUS_TAXON = "Q763978"
Q_CANDIDATUS = "Q857968"
Q_INC_SED = "Q21397852"
Q_WASTEBASKET = "Q962530"
Q_UNDESCRIBED = "Q7883954"
Q_SUBTYPE = "Q19862317"
Q_STRAIN = "Q855769"
Q_EXTINCT_TAXON = "Q98961713"
Q_TAX_HYP = "Q124477390"
TAXA_SET = frozenset(
    [
        Q_TAXON,
        Q_FOSSIL_TAXON,
        Q_MONOTYPIC_FOSSIL_TAXON,
        Q_CLADE,
        Q_MONOTYPIC_TAXON,
        Q_INC_SED,
        Q_SPECIES,
        Q_SUBSPECIES,
        Q_MONOPHYLY,
        Q_LAZARUS_TAXON,
        Q_CANDIDATUS,
        Q_WASTEBASKET,
        Q_UNDESCRIBED,
        Q_SUBTYPE,
        Q_STRAIN,
        Q_EXTINCT_TAXON,
        Q_TAX_HYP,
    ]
)


def wikid_ent_is_taxon(item: WikidataItem) -> bool:
    """Return True if the Wikidata Item has occupation politician."""
    claim_group = item.get_claim_group(P_IS_INSTANCE_OF)
    instance_qids = set(
        [
            claim.mainsnak.datavalue.value["id"]
            for claim in claim_group
            if claim.mainsnak.snaktype == "value"
        ]
    )
    ret = not TAXA_SET.isdisjoint(instance_qids)
    # _LOG.debug(f"enitity {item.entity_id} P_IS_INSTANCE_OF {instance_qids} IS_TAXON={ret}")
    return ret


def _parse_taxonomy_file(taxonomy_fp):
    lp = tax_wikidata_parser
    id_2_taxon = {}
    with open(taxonomy_fp, "r") as inp:
        lit = iter(inp)
        fl = next(lit)
        assert fl == TAXWIKIDATA_HEADER
        for n, line in enumerate(lit):
            ls = line.strip()
            if not ls:
                continue
            if n % 10000 == 0 and n > 0:
                _LOG.debug(' read taxon {:<7} from "{}" ...'.format(n, taxonomy_fp))
            obj = Taxon(line, line_parser=lp)
            if obj.id in id_2_taxon:
                _LOG.warning(f"Duplicate taxon ID: {obj.id}")
                if obj.__dict__ != id_2_taxon[obj.id].__dict__:
                    m = f"Duplicate taxon ID: {obj.id} with differing content."
                    raise RuntimeError(m)
            id_2_taxon[obj.id] = obj
    return id_2_taxon


EMIT_ID_FOR_PROP = frozenset(
    [
        # WARN
        "P13177",  # "homonymous taxon" - symmetrical
        # FLAG hybrid
        "P1531",  # "hybrid of",
    ]
)
OBJ_IS_SYN = frozenset(
    [
        "P1420",
        "P1403",  # "original combination", - complement of basionym or protonym
    ]
)
ENT_IS_SYN = frozenset(
    [
        "P12763",  # "taxon synonym of" - found in jr synonym listing valid name
        "P12764",  # "replaced synonym of", - found in jr synonym listing valid name
        "P694",  # "replaced synonym (for nom. nov.)",
    ]
)
PRED_IS_IGNORED = frozenset(
    [
        "P12765",  # "protonym of", -  found in jr synonym listing valid name
        "P12766",  # "basionym of", -  found in jr synonym listing valid name
        "P13177",  # "homonymous taxon" - symmetrical
        # IGNORE for now
        "P13478",  # "nomenclatural type of",
        "P427",  # "taxonomic type"
        # IGNORE for now
        "P5304",  # "type locality (biology)",
        "P1137",  # "fossil found in this unit" - erroneous subject should be strat. layer
    ]
)


def read_additional_props(id_2_taxon, additional_props_fp):
    synonym_set = set()
    with open(additional_props_fp, "r") as inp:
        lit = iter(inp)
        fl = next(lit)
        assert fl == "Entity\tPredicate\tObject\n"
        for line in lit:
            ls = [i.strip() for i in line.split("\t")]
            e_id, pred, obj_id = ls
            try:
                taxon = id_2_taxon[e_id]
            except:
                _LOG.warning(f"Taxon {e_id} not among taxa!")
                continue
            if pred in PRED_IS_IGNORED:
                continue
            if pred in ENT_IS_SYN:
                synonym_set.add(e_id)
                taxon.add_synonym_id(obj_id)
            elif pred in OBJ_IS_SYN:
                if obj_id not in id_2_taxon:
                    _LOG.warning(f"synonym object {obj_id} for {e_id} not among taxa!")
                else:
                    synonym_set.add(obj_id)
                    syn_taxon = id_2_taxon[obj_id]
                    syn_taxon.add_synonym_id(e_id)
            elif pred == "P2093":
                taxon.author_str = obj_id
            elif pred == "P1531":
                taxon.flag_as_hybrid()
            else:
                raise RuntimeError(f"Unexpected predicate {pred}")
    return synonym_set


def parse_wikidata(taxonomy_fp, additional_props_fp):
    id_2_taxon = _parse_taxonomy_file(taxonomy_fp)
    synonyms = read_additional_props(id_2_taxon, additional_props_fp)
    return id_2_taxon, synonyms
