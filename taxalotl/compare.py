#!/usr/bin/env python
from __future__ import print_function
from peyotl import get_logger
import sys
import os
from taxalotl.tax_partition import (get_taxon_partition, INP_TAXONOMY_DIRNAME, MISC_DIRNAME, OUTP_TAXONOMY_DIRNAME)
from .util import OutFile, OutDir
_LOG = get_logger(__name__)


def get_frag_from_dir(taxalotl_conf, tax_dir):
    res = taxalotl_conf.get_terminalized_res_by_id("ott")
    pd = res.partitioned_filepath
    assert tax_dir.startswith(pd)
    f = tax_dir[len(pd):]
    while f.startswith('/'):
        f = f[1:]
    return f


def compare_taxonomies_in_dir(taxalotl_conf, tax_dir):
    fragment = get_frag_from_dir(taxalotl_conf, tax_dir)
    _LOG.info("fragment = {}".format(fragment))
    tax_id_set = set()
    non_misc_dir = os.path.join(tax_dir, INP_TAXONOMY_DIRNAME)
    misc_dir = os.path.join(tax_dir, MISC_DIRNAME, INP_TAXONOMY_DIRNAME)
    for sd in [misc_dir, non_misc_dir]:
        if os.path.exists(sd):
            tax_id_set.update(os.listdir(sd))
    out = sys.stdout
    out_dir = os.path.join(tax_dir, OUTP_TAXONOMY_DIRNAME)
    with OutDir(out_dir):
        for res_id in tax_id_set:
            res = taxalotl_conf.get_resource_by_id(res_id)
            tp = get_taxon_partition(res, fragment)
            tp.read_inputs_for_read_only()
            tf = tp.get_taxa_as_forest()
            fn = os.path.split(fragment)[-1]
            fp = os.path.join(out_dir, '{}-for-{}.txt'.format(res_id, fn))
            with OutFile(fp) as outstream:
                out.write("writing taxonomy for {} according to {} to \"{}\"\n".format(fragment, res_id, fp))
                tf.write_indented(outstream)
            semantics_dir = os.path.join(out_dir, res_id)
            with OutDir(semantics_dir):
                res.semanticize(fragment, semantics_dir, tax_part=tp, taxon_forest=tf)
