#!/usr/bin/env python
import sys
import os


def _do_name_grep_taxonomy(fp, name_pat, outstream=sys.stdout):
    return _do_grep_of_col(fp, name_pat, 2, outstream=outstream)


def _do_name_grep_synonyms(fp, name_pat, outstream=sys.stdout):
    return _do_grep_of_col(fp, name_pat, 0, outstream=outstream)


def _do_tax_id_grep_taxonomy(fp, tax_id_pat, outstream=sys.stdout):
    return _do_grep_of_col(fp, tax_id_pat, 0, outstream=outstream)


def _do_tax_id_grep_synonyms(fp, tax_id_pat, outstream=sys.stdout):
    return _do_grep_of_col(fp, tax_id_pat, 1, outstream=outstream)


def _do_grep_of_col(fp, name_pat, col_idx, outstream=sys.stdout):
    if not os.path.exists(fp):
        return []
    matches = []
    with open(fp, "r") as inp:
        li = iter(inp)
        next(li)
        for line in li:
            ls = line.split("\t|\t")
            name_str = ls[col_idx]
            if name_pat.match(name_str):
                matches.append(line[:-1])
    if not matches:
        return []
    tp = fp
    if "partitioned" in tp:
        tp = tp.split("partitioned")[-1]
    while tp.startswith("/"):
        tp = tp[1:]
    while tp.endswith("/"):
        tp = tp[:-1]
    pref = f"{tp}:"
    ret = []
    for line in matches:
        ret.append((fp, line))
        if outstream:
            outstream.write(f"{pref} {line}\n")
    return ret


def grep_name_in_single_res(
    rw, name_pat, target, outstream=sys.stdout, root_taxon=None
):
    return _generic_grep_in_one_res(
        rw,
        name_pat,
        tax_fn=_do_name_grep_taxonomy,
        syn_fn=_do_name_grep_synonyms,
        target=target,
        outstream=outstream,
        root_taxon=root_taxon,
    )


def grep_tax_id_in_single_res(
    rw, tax_id_pat, target, outstream=sys.stdout, root_taxon=None
):
    return _generic_grep_in_one_res(
        rw,
        tax_id_pat,
        tax_fn=_do_tax_id_grep_taxonomy,
        syn_fn=_do_tax_id_grep_synonyms,
        target=target,
        outstream=outstream,
        root_taxon=root_taxon,
    )


_SEARCH_TAX_SET = frozenset(["both", "taxa"])
_SEARCH_SYN_SET = frozenset(["both", "synonyms"])


def grep_name_in_single_res_taxa(rw, name_pat, outstream=sys.stdout, root_taxon=None):
    return grep_name_in_single_res(
        rw, name_pat, target="taxa", outstream=outstream, root_taxon=root_taxon
    )


def grep_name_in_single_res_syn(rw, name_pat, outstream=sys.stdout, root_taxon=None):
    return grep_name_in_single_res(
        rw, name_pat, target="synonyms", outstream=outstream, root_taxon=root_taxon
    )


def grep_tax_id_in_single_res_taxa(
    rw, tax_id_pat, outstream=sys.stdout, root_taxon=None
):
    return grep_tax_id_in_single_res(
        rw, tax_id_pat, target="taxa", outstream=outstream, root_taxon=root_taxon
    )


def grep_tax_id_in_single_res_syn(
    rw, tax_id_pat, outstream=sys.stdout, root_taxon=None
):
    return grep_tax_id_in_single_res(
        rw, tax_id_pat, target="synonyms", outstream=outstream, root_taxon=root_taxon
    )


def _generic_grep_in_one_res(
    rw, pat, tax_fn, syn_fn, target="both", outstream=sys.stdout, root_taxon=None
):
    search_tax = target.lower() in _SEARCH_TAX_SET
    search_syn = target.lower() in _SEARCH_SYN_SET
    if not (search_tax or search_syn):
        msg = f"target should be both, taxa, or synonyms. Got {target}"
        raise ValueError(msg)
    r = []
    if rw.has_been_partitioned:
        if search_tax:
            tfp = rw.get_part_taxa_filepaths(below=root_taxon)
            for fn in tfp:
                r.extend(tax_fn(fn, pat, outstream=outstream))
        if search_syn:
            sfp = rw.get_part_syn_filepaths(below=root_taxon)
            for fn in sfp:
                r.extend(syn_fn(fn, pat, outstream=outstream))
        return r
    if not rw.has_been_normalized:
        raise RuntimeError(
            f"{rid} needs to be normalized or partitioned to work with grep"
        )
    if search_tax:
        fn = os.path.join(rw.normalized_filedir, "taxonomy.tsv")
        r = tax_fn(fn, pat, outstream=outstream)
    if search_syn:
        fn = os.path.join(rw.normalized_filedir, "synonyms.tsv")
        r.extend(syn_fn(fn, pat, outstream=outstream))
    return r
