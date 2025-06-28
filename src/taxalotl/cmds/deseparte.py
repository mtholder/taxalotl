#!/usr/bin/env python
from __future__ import print_function
import logging
import sys
import os

from ..tax_partition import (
    TAXONOMY_FN,
    SYNONYMS_FN,
    ROOTS_FILENAME,
    INP_TAXONOMY_DIRNAME,
    MISC_DIRNAME,
    OUTP_TAXONOMY_DIRNAME,
)
from ..util import OutFile, OutDir, get_frag_from_dir

_LOG = logging.getLogger(__name__)


def deseparate_taxonomies_in_dir(taxalotl_conf, tax_dir):
    fragment = get_frag_from_dir(taxalotl_conf, tax_dir)
    _LOG.info("fragment = {}".format(fragment))
    tax_id_set = set()
    misc_dir = os.path.join(tax_dir, MISC_DIRNAME, INP_TAXONOMY_DIRNAME)
    par_dir = os.path.split(tax_dir)[0]
    if os.path.isdir(misc_dir):
        m = "The presence of a directory at {} indicates that deseparate is not a safe command."
        raise RuntimeError(m.format(misc_dir))
    par_inp = os.path.join(par_dir, INP_TAXONOMY_DIRNAME)
    par_misc = os.path.join(par_dir, MISC_DIRNAME, INP_TAXONOMY_DIRNAME)
    non_misc_dir = os.path.join(tax_dir, INP_TAXONOMY_DIRNAME)
    if os.path.exists(non_misc_dir):
        tax_id_set.update(os.listdir(non_misc_dir))
    out = sys.stdout
    for res_id in tax_id_set:
        inp_dir = os.path.join(non_misc_dir, res_id)
        tax_fp = os.path.join(inp_dir, TAXONOMY_FN)
        syn_fp = os.path.join(inp_dir, SYNONYMS_FN)
        roots_fp = os.path.join(inp_dir, ROOTS_FILENAME)
        if not any([os.path.isfile(i) for i in [tax_fp, syn_fp, roots_fp]]):
            continue
        tax_dest, syn_dest = None, None
        m = "The absence of a destination for {} in {} made the deseprate command fail for {}"
        if os.path.isfile(tax_fp):
            tax_dest = _find_destination_fp(res_id, TAXONOMY_FN, par_inp, par_misc)
            if tax_dest is None:
                _LOG.warning(m.format(tax_fp, par_inp, res_id))
                continue
        if os.path.isfile(syn_fp):
            syn_dest = _find_destination_fp(res_id, TAXONOMY_FN, par_inp, par_misc)
            if syn_dest is None:
                _LOG.warning(m.format(syn_fp, par_inp, res_id))
                continue
        if tax_dest:
            _move_all_but_first_line(tax_fp, tax_dest)
        if syn_dest:
            _move_all_but_first_line(syn_fp, syn_dest)
        for i in [tax_fp, syn_fp, roots_fp]:
            if os.path.isfile(i):
                _LOG.debug("Removing {}".format(i))
                os.remove(i)
        try:
            _LOG.debug("Removing {}".format(inp_dir))
            os.rmdir(inp_dir)
        except:
            _LOG.exception('Could not remove dir "{}"'.format(inp_dir))


def _find_destination_fp(res_id, fn, dir1, dir2):
    for d in [dir1, dir2]:
        fp = os.path.join(d, res_id, fn)
        if os.path.isfile(fp):
            return fp
    return None


def _move_all_but_first_line(inp_fp, dest_fp):
    _LOG.debug('Copying content from "{}"  to "{}"'.format(inp_fp, dest_fp))
    seen_lines = open(dest_fp, "r", encoding="utf-8").readlines()
    seen_lines = set(seen_lines)
    with open(inp_fp, "r", encoding="utf-8") as inp:
        with OutFile(dest_fp, mode="a") as outp:
            for n, line in enumerate(inp):
                if n == 0:
                    assert line.startswith("uid") or line.startswith("name\t|\tuid")
                else:
                    if line not in seen_lines:
                        outp.write(line)
