#!/usr/bin/env python

# Based on example script `basic_json_dump.py
# from qwikidata
#     https://qwikidata.readthedocs.io/en/stable/readme.html
import sys
from qwikidata.entity import WikidataItem
from qwikidata.json_dump import WikidataJsonDump

P_IS_INSTANCE_OF = "P31"
Q_TAXON = "Q16521"
Q_FOSSIL_TAXON = "Q23038290"

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
                continue
        else:
            sys.stderr.write(f"Skipping non taxon {entity.entity_id} ...\n")
            continue
        taxa.append(entity)
    sys.stderr.write(f"{len(taxa)} taxa found\n")

if __name__ == "__main__":
    sys.exit(main(inp_fp=sys.argv[1]))