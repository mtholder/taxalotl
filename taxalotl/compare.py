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
    graph_by_res_id = {}
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
                graph = res.semanticize(fragment, semantics_dir, tax_part=tp, taxon_forest=tf)
                graph_by_res_id[res_id] = (res, graph)
    ott_res = taxalotl_conf.get_terminalized_res_by_id("ott", None)
    ott_graph = graph_by_res_id[ott_res.id][1]
    ott_vstc = ott_graph.valid_specimen_based_taxa
    ott_vn2tc = ott_graph.valid_name_to_taxon_concept_map
    for res_id, res_graph_pair in graph_by_res_id.items():
        if res_id == ott_res.id:
            continue
        res, ref_graph = res_graph_pair
        out.write('comparing {} to {}\n'.format(ott_res.id, res_id))
        ref_vstc = ref_graph.valid_specimen_based_taxa
        out.write('{} vs {} valid specimen-based names\n'.format(len(ott_vstc), len(ref_vstc)))
        ref_vn2tc = ref_graph.valid_name_to_taxon_concept_map
        just_ott, both, just_ref =[], [], []
        for ott_name, tax_con in ott_vn2tc.items():
            if not tax_con.is_specimen_based:
                continue
            t = both if ott_name in ref_vn2tc else just_ott
            t.append(ott_name)
        for ott_name in ref_vn2tc.keys():
            if ott_name not in ott_vn2tc:
                just_ref.append(ott_name)
        just_ott.sort()
        just_ref.sort()
        _write_just_in_list(out, just_ott, ott_res, ott_vn2tc)
        _write_just_in_list(out, just_ref, res, ref_vn2tc)


def _write_just_in_list(out, just_in, res, obj_lookup):
    out.write('{} only in  {} :\n'.format(len(just_in), res.id))
    for n, i in enumerate(just_in):
        out.write('{} "{}" : '.format(1 + n, i))
        obj = obj_lookup[i]
        obj.explain(out)
        out.write('\n')
