# Tweaked example script from qwikidata
import sys
import time
from qwikidata.entity import WikidataItem
from qwikidata.json_dump import WikidataJsonDump
from qwikidata.utils import dump_entities_to_json
from subprocess import call

P_IS_INSTANCE_OF = "P31"
Q_TAXON = "Q16521"
Q_FOSSIL_TAXON = "Q23038290"
Q_MONOTYPIC_FOSSIL_TAXON = "Q47487597"
Q_CLADE = "Q713623"
Q_MONOTYPIC_TAXON = "Q310890"
Q_SPECIES = "Q15458879"
Q_SUBSPECIES = "Q68947"
Q_MONOPHYLY = "Q210958" # Sigh...
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
DEBUG = "-v" in sys.argv


def debug(line):
    if DEBUG:
        sys.stderr.write(f"{line}\n")


def is_taxon(item: WikidataItem) -> bool:
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
    debug(f"enitity {item.entity_id} P_IS_INSTANCE_OF {instance_qids} IS_TAXON={ret}")
    return ret


def dump_block_of_taxa(taxa, out_fp):
    sys.stderr.write(f"Dumping {len(taxa)} to {out_fp} ...\n")
    dump_entities_to_json(taxa, out_fp)
    retcode = call(["bzip2", out_fp])
    if retcode != 0:
        raise RuntimeError(f"bzip2 of {out_fp} failed.\n")


def main(inp_fp, out_fp):
    # create an instance of WikidataJsonDump
    wjd_dump_path = sys.argv[1]
    wjd = WikidataJsonDump(wjd_dump_path)

    # create an iterable of WikidataItem representing politicians
    taxa = []
    t1 = time.time()
    prev_dumped = 0
    for ii, entity_dict in enumerate(wjd):
        if entity_dict["type"] == "item":
            entity = WikidataItem(entity_dict)
            if is_taxon(entity):
                taxa.append(entity)
        nt = len(taxa)
        if ii % 1000 == 0:
            t2 = time.time()
            dt = t2 - t1
            rate = ii / dt
            tt = nt + prev_dumped
            sys.stderr.write(
                f"found {tt} taxa among {ii} entities [entities/s: {rate:.2f}]\n"
            )
        if nt == 1000:
            prev_dumped += nt
            dump_block_of_taxa(taxa, f"block-{prev_dumped}-{out_fp}")
            taxa = []
    nt = len(taxa)
    prev_dumped += nt
    sys.stderr.write(f"found {prev_dumped} taxa among {ii} entities\n")
    if taxa:
        dump_block_of_taxa(taxa, f"block-{prev_dumped}-{out_fp}")
    sys.stderr.write(f"Done!\n")


if __name__ == "__main__":
    sys.exit(main(inp_fp=sys.argv[1], out_fp="filtered_entities.json"))
