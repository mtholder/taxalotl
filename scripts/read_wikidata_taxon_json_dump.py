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
P_TAXON_RANK = "P105"

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
    "Q2372065": "nomen correctum",
    "Q123561351": "ineditus",
    "Q17134993": "validly published name",
    "Q17487588": "unavailable combination",
    "Q7882332": "unavailable name",
    "Q15709300": "nomen protectum",
    "Q28549151": "preoccupied name",
    "Q1040689": "synonym",
    "Q42106": "synonym", # should be corrected to Q1040689
    "Q3766304": "species inquirenda",
    "Q17276484": "later homonym",
    "Q2491016": "orthographical variant",
    "Q18912752": "disputed",
    "Q17086880": "valid",
    "Q14594740": "recombination",
    "Q67943587": "genus neutrum conservandum",
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


DIE_ON = frozenset([
  "P1135", 
  "P678", 
])

EMIT_ID_FOR_PROP = frozenset([
    "P1137",  # "fossil found in this unit"
    "P12763", # "taxon synonym of",
    "P12764", # "replaced synonym of",
    "P12765", # "protonym of",
    "P12766", # "basionym of",
    "P13177", # "homonymous taxon"
    "P13478", # "nomenclatural type of",
    "P1531",  # "hybrid of",
    "P1403",  # "original combination",
    "P1420",  # "taxon synonym"
    "P427",   # "taxonomic type"
    "P5304",  # "type locality (biology)",
    "P694",   # "replaced synonym (for nom. nov.)",  
])

EMIT_STR_FOR_PROP = frozenset([
    "P2093", # "author name string",
])

TO_STD_ERROR = frozenset([
  "P2241",  
  "P7452",  ])

def extra_log(line):
    out = sys.stderr
    out.write(f"{line}\n")

def warn(line):
    sys.stderr.write(f"WARNING: {line}\n")

def emit_ids_for_prop(eid, pid, claims_group):
    for claim in claims_group:
        dv = claim.mainsnak.datavalue
        if dv is None:
            warn(f"{eid} {pid} has no datavalue")
            continue
        try:
            id_val = dv.value["id"]
        except:
            raise RuntimeError(f"error extracting ID for {pid} for entity {eid}")
        extra_log(f"{eid}\t{pid}\t{id_val}")

def emit_str_for_prop(eid, pid, claims_group):
    for claim in claims_group:
        dv = claim.mainsnak.datavalue
        if dv is None:
            warn(f"{eid} {pid} has no str datavalue")
            continue
        try:
            s_val = str(dv.value)
        except:
            raise RuntimeError(f"error extracting str for {pid} for entity {eid}")
        extra_log(f"{eid}\t{pid}\t{s_val}")

def format_src(ext_id):
    if not ext_id:
        return ""
    full_list = []
    for k, vlist in ext_id.items():
        plist = [f"{k}:{v}" for v in vlist]
        full_list.append(",".join(plist))
    return ",".join(full_list)

ENT_TO_RANK_NAME = {
    "Q113015256" : "ichnospecies",
    "Q1153785" : "subphylum",
    "Q1306176" : "nothospecies",
    "Q13198444" : "subseries",
    "Q14817220" : "supertribe",
    "Q164280" : "subfamily",
    "Q19858692" : "superkingdom",
    "Q2007442" : "infraclass",
    "Q21061204" : "subterclass",
    "Q21074316" : "hyporder",
    "Q2111790" : "superphylum",
    "Q2136103" : "superfamily",
    "Q227936" : "tribe",
    "Q2361851" : "infraphylum",
    "Q2455704" : "subfamily",
    "Q2752679" : "subkingdom",
    "Q279749" : "form",
    "Q2889003" : "infraorder",
    "Q2981883" : "cohort",
    "Q3025161" : "series", # botany
    "Q3181348" : "section",
    "Q3238261" : "subgenus",
    "Q334460" : "division",
    "Q34740" : "genus",
    "Q3504061" : "superclass",
    "Q35409" : "family",
    "Q36602" : "order",
    "Q36732" : "kingdom",
    "Q37517" : "class",
    "Q3825509" : "forma specialis",
    "Q38348" : "phylum",
    "Q3965313" : "subtribe",
    "Q4150646" : "cultivar group",
    "Q4886" : "cultivar",
    "Q5867051" : "subclass",
    "Q5867959" : "suborder",
    "Q5868144" : "superorder",
    "Q5998839" : "subsection", # botany
    "Q630771" : "subvariety",
    "Q6311258" : "parvorder",
    "Q6541077" : "subcohort",
    "Q68947" : "subspecies",
    "Q713623" : "clade",
    "Q7432" : "species",
    "Q7506274" : "mirorder",
    "Q767728" : "variety",
    "Q10861375": "subsection", # zoology
    "Q3491997": "subdivision",
    "Q21061732": "series", # zoology
    "Q6045742": "nothogenus",
    "Q6054535": "infralegion",
    "Q6054637": "sublegion",
    "Q7504331": "legion",
    "Q6054795": "superlegion",
    "Q10861426": "section",
    "Q3150876": "infrakingdom",
    "Q23760204": "superdivision",
    "Q59779184": "nothovariety",
    "Q146481": "domain",
    "Q59772023": "nothosubspecies",
    "Q6462265": "grandorder",
    "Q5469884": "form",
    "Q112082101": "ichnogenus",
    "Q1972414": "pathovar",
    "Q62075839": "realm",
    "Q7574964": "species group",
}

def to_rank_name(rank_q):
    return ENT_TO_RANK_NAME[rank_q]

def process_taxon(taxon):
    eid = taxon.entity_id
    claims = taxon.get_truthy_claim_groups()
    # sys.exit(f"Claims =\n  {repr(claims)}\n")
    parent, tax_name, tax_aut, tax_year, rank = "", "", "", "", ""
    ext_id = {}
    nom_status = None
    for pid, wcg in claims.items():
        if pid in EMIT_ID_FOR_PROP:
            emit_ids_for_prop(eid, pid, wcg)
            continue
        elif pid in EMIT_STR_FOR_PROP:
            emit_str_for_prop(eid, pid, wcg)
            continue
        if pid in DIE_ON:
            raise RuntimeError(f"prop {pid} in DIE_ON list")
        if pid in TO_STD_ERROR:
            raise RuntimeError(f"prop {pid} in TO_STD_ERROR list")
        if pid == P_TAXON_RANK:
            dv_list = [c.mainsnak.datavalue for c in wcg if c.mainsnak.datavalue is not None]
            ilist = [to_rank_name(dv.value["id"]) for dv in dv_list]
            if ilist:
                rank = ",".join(ilist)
            continue
        for claim in wcg:
            dv = claim.mainsnak.datavalue
            # pid = claim.property_id
            if pid == P_PARENT:
                try:
                    parent = dv.value["id"]
                except AttributeError:
                    warn(f"NO PARENT TAXON for {eid} !")
            elif pid == P_TAXON_NAME:
                try:
                    tax_name = str(dv.value)
                except AttributeError:
                    warn(f"No tax name for {eid} !")
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
                        m = f"unknown qual {qid} in tax_name for {eid}"
                        warn(m)
            elif pid in EXT_ID_TO_PREF:
                if dv is not None:
                    pref = EXT_ID_TO_PREF[pid]
                    if pref in ext_id:
                        ext_id[pref].append(str(dv.value))
                    else:
                        ext_id[pref] = [str(dv.value)]

    if not nom_status:
        nom_status = ""
    elif isinstance(nom_status, list):
        nom_status = ",".join(nom_status)
    src_str = format_src(ext_id)

    aut_nf = ",".join(tax_aut) if isinstance(tax_aut, list) else tax_aut
    aut_yf = ",".join(tax_year) if isinstance(tax_year, list) else tax_year
    el_list = [eid, 
               parent,
               str(tax_name),
               rank,
               aut_nf,
               aut_yf,
               nom_status,
               src_str]
    sep = "\t|\t"
    print(sep.join(el_list))


def wr_process_taxon(taxon):
    try:
        return process_taxon(taxon)
    except:
        df = json.dumps(taxon._entity_dict, sort_keys=True, indent=2)
        warn(f"ERROR on:\n{df}")
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
                warn(f"Skipping non taxon {entity.entity_id} at idx={idx}...")
            else:
                taxa.append(entity)
        else:
            warn(f"Skipping non taxon {entity.entity_id} ...")
        
    # warn(f"{len(taxa)} taxa found\n")
    for taxon in taxa:
        wr_process_taxon(taxon)

if __name__ == "__main__":
    sys.exit(main(inp_fp=sys.argv[1]))