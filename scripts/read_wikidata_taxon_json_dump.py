#!/usr/bin/env python

# Based on example script `basic_json_dump.py
# from qwikidata
#     https://qwikidata.readthedocs.io/en/stable/readme.html
import json
import sys
from qwikidata.entity import WikidataItem
from qwikidata.json_dump import WikidataJsonDump

P_IS_INSTANCE_OF = "P31"
Q_TAXON = "Q16521"
Q_FOSSIL_TAXON = "Q23038290"

P_PARENT = "P171"
P_TAXON_NAME = "P225"
P_TAXON_AUTHOR = "P405"
P_TAXON_YEAR = "P574"

P_OTT = "P9157"
P_GBIF = "P846"
P_COL = "P10585"
P_NCBI = "P685"
P_IRMNG = "P5055"

P_NOM_STAT = "P1135"

P_FOSSIL_FOUND = "P1137"
P_TAX_SYN_OF = "P12763"
P_REPLACED_SYN_OF = "P12764"
P_PROTONYM_OF = "P12765"
P_BASIONYM_OF = "P12766"
P_HOMONYM = "P13177"
P_NOM_TYPE_OF = "P13478"
P_TAX_SYN = "P1420"
P_HYBRID_OF = "P1531"
P_AUTH_STR = "P2093"
P_INC_SED = "P678"
P_REP_SYN_NOM_NOV = "P694"
P_TYPE_LOC = "P5304"
  
def is_taxon(item: WikidataItem, truthy: bool = True) -> bool:
    """Return True if the Wikidata Item is a taxon or fossil taxon instance"""
    if truthy:
        claim_group = item.get_truthy_claim_group(P_IS_INSTANCE_OF)
    else:
        claim_group = item.get_claim_group(P_IS_INSTANCE_OF)

    is_a_qids = [
        claim.mainsnak.datavalue.value["id"]
        for claim in claim_group
        if claim.mainsnak.snaktype == "value"
    ]
    return (Q_TAXON in is_a_qids) or (Q_FOSSIL_TAXON in is_a_qids)

def get_year(qual_prop):
    qdv = qual_prop.snak.datavalue
    try:
        raw = qdv.value["time"]
    except:
        return None
    cm = qdv.value["calendarmodel"]
    assert(cm == "http://www.wikidata.org/entity/Q1985727")
    assert(raw[0] == "+")
    return raw[1:5]

KNOWN_TAX_NAME_QUAL = frozenset([
    P_TAXON_AUTHOR, 
    P_TAXON_YEAR,])

SKIPPABLE_TN_QUAL = frozenset([
    "P3831", # "object of statement has role"
    "P697",  # "ex taxon author"
    "P2433", # "gender of sci name"
    "P1353", # "original spelling"
    "P138",  # "named after"
    ])

NOM_STAT_MAP = {
    "Q1093954": "nomen illegitimum",
    "Q30349290": "nomen invalidum",
    "Q15149791": "nomen utique rejiciendum",
    "Q941227": "nomen conservandum",
    "Q922448": "nomen dubium",
    "Q29995613": "superfluous name",
    "Q51165454": "typus conservandus",
    "Q68455470": "protected name",
    "Q28597224": "needing care",
    "Q51165361": "orthographia conservanda",
    "Q23038368": "nomen rejiciendum propositum",
    "Q17276482": "nomen rejiciendum",
    "Q749462": "replacement name",
    "Q15708833": "nomen conservandum propositum",
    "Q844326": "nomen nudum",
    "Q122735442": "nomen praeoccupatum",
    "Q130297303": "not validly published name",
    "Q15708837": "nomen utique rejiciendum propositum",
    "Q901847": "nomen oblitum",
    "Q70573632": "nomen approbatum",
    }

EXT_ID_TO_PREF = {
    P_OTT : "ott",
    P_GBIF : "gbif",
    P_COL : "col",
    P_NCBI : "ncbi",
    P_IRMNG : "irmng",    
}

def _get_tax_aut(qlist):
    assert(len(qlist) > 0)
    ret = []
    for au in qlist:
        try:
            ret.append(au.snak.datavalue.value["id"])
        except AttributeError:
            pass
    return ret

def get_nom_status(qlist):
    ilist = [q.snak.datavalue.value["id"] for q in qlist]
    slist = [NOM_STAT_MAP[i] for i in ilist]
    if len(slist) == 1:
        return slist[0]
    return slist

def get_tn_year(qlist):
    ys = set([get_year(q) for q in qlist])
    yl = [i for i in ys if i is not None]
    if len(yl) == 1:
        return yl[0]
    return yl

def process_taxon(taxon):
    eid = taxon.entity_id
    claims = taxon.get_truthy_claim_groups()
    # sys.exit(f"Claims =\n  {repr(claims)}\n")
    parent, tax_name, tax_aut, tax_year = "", "", "", ""
    ext_id = {}
    nom_status = None
    for pid, wcg in claims.items():
        for claim in wcg:
            dv = claim.mainsnak.datavalue
            # pid = claim.property_id
            if pid == P_PARENT:
                try:
                    parent = dv.value["id"]
                except AttributeError:
                    sys.stderr.write(f"NO PARENT TAXON for {eid}!\n")
            elif pid == P_TAXON_NAME:
                tax_name = str(dv)
                for qid, qlist in claim.qualifiers.items():
                    if qid in SKIPPABLE_TN_QUAL:
                        continue
                    if qid == P_TAXON_AUTHOR:
                        tax_aut = _get_tax_aut(qlist)
                    elif qid == P_TAXON_YEAR:
                        tax_year = get_tn_year(qlist)
                    elif qid == P_NOM_STAT:
                        nom_status = get_nom_status(qlist)
                    else:
                        m = f"WARN: unkown qual {qid} in tax_name for {eid}\n"
                        sys.stderr.write(m)
            elif pid in EXT_ID_TO_PREF:
                pref = EXT_ID_TO_PREF[pid]
                if pref in ext_id:
                    ext_id[pref].append(str(dv))
                else:
                    ext_id[pref] = [str(dv)] 


def wr_process_taxon(taxon):
    try:
        return process_taxon(taxon)
    except:
        df = json.dumps(taxon._entity_dict, sort_keys=True, indent=2)
        sys.stderr.write(f"ERROR on:\n{df}\n")
        raise

def main(inp_fp):
    # create an instance of WikidataJsonDump
    wjd_dump_path = sys.argv[1]
    wjd = WikidataJsonDump(wjd_dump_path)

    # create an iterable of WikidataItem representing politicians
    taxa = []
    for idx, entity_dict in enumerate(wjd):
        if entity_dict["type"] == "item":
            entity = WikidataItem(entity_dict)
            if not is_taxon(entity):
                sys.stderr.write(f"Skipping non taxon {entity.entity_id} at idx={idx}...\n")
            else:
                taxa.append(entity)
        else:
            sys.stderr.write(f"Skipping non taxon {entity.entity_id} ...\n")
        
    # sys.stderr.write(f"{len(taxa)} taxa found\n")
    for taxon in taxa:
        wr_process_taxon(taxon)

if __name__ == "__main__":
    sys.exit(main(inp_fp=sys.argv[1]))