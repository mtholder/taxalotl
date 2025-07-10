# Tweaked example script from qwikidata
import sys
import time
from qwikidata.entity import WikidataItem
from qwikidata.json_dump import WikidataJsonDump
from qwikidata.utils import dump_entities_to_json
from subprocess import call
from taxalotl.wikidata import wikid_ent_is_taxon


def dump_block_of_taxa(taxa, out_fp):
    sys.stderr.write(f"Dumping {len(taxa)} to {out_fp} ...\n")
    dump_entities_to_json(taxa, out_fp)
    retcode = call(["bzip2", out_fp])
    if retcode != 0:
        raise RuntimeError(f"bzip2 of {out_fp} failed.\n")


def main(inp_fp, out_fp):
    # create an instance of WikidataJsonDump
    wjd = WikidataJsonDump(inp_fp)

    taxa = []
    t1 = time.time()
    prev_dumped = 0
    ii = 0
    for ii, entity_dict in enumerate(wjd):
        if entity_dict["type"] == "item":
            entity = WikidataItem(entity_dict)
            if wikid_ent_is_taxon(entity):
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
