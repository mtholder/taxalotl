#!/usr/bin/env python
import sys


def main(taxonomy_fp, par_id_fp):
    with open(par_id_fp, "r") as inp:
        id_list = [i.strip() for i in inp.readlines() if i.strip()]
    par_ids = frozenset(id_list)
    with open(taxonomy_fp, "r") as inp:
        for line in inp.readlines():
            ls = line.strip().split("\t")
            assert ls[1] == "|"
            par_id = ls[2]
            if par_id in par_ids:
                print(line[:-1])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
