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

JUST_COF = True

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
            if JUST_COF and not (res_id.startswith('cof') or res_id.startswith('ott')):
                continue
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
    for res_id, res_graph_pair in graph_by_res_id.items():
        if res_id == ott_res.id:
            continue
        if JUST_COF and not res_id.startswith('cof'):
            continue
        res, ref_graph = res_graph_pair
        _write_report(out, ott_res.id, ott_graph, res_id, ref_graph)

def _write_report(out, ott_id, ott_graph, res_id, ref_graph):
    ott_vstc = ott_graph.valid_specimen_based_taxa
    ott_vn2tc = ott_graph.valid_name_to_taxon_concept_map
    out.write('comparing {} to {}\n'.format(ott_id, res_id))
    ref_vstc = ref_graph.valid_specimen_based_taxa
    ref_vn2tc = ref_graph.valid_name_to_taxon_concept_map
    just_ott, both, just_ref = [], [], []
    for ott_name, tax_con in ott_vn2tc.items():
        if not tax_con.is_specimen_based:
            continue
        ref_tc = ref_vn2tc.get(ott_name)
        if ref_tc:
            both.append((ott_name, ott_name, tax_con, ref_tc))
        else:
            just_ott.append(ott_name)
    for ott_name in ref_vn2tc.keys():
        if ott_name not in ott_vn2tc:
            just_ref.append(ott_name)
    just_ott.sort()
    just_ref.sort()
    out.write('{} only in  {} :\n'.format(len(just_ott), ott_id))
    dnm = _write_report_info_unmatched(out, just_ott, ott_vn2tc, ref_graph)
    diff_name_matched = dnm
    out.write('{} only in  {} :\n'.format(len(just_ref), res_id))
    dnm = _write_report_info_unmatched(out, just_ref, ref_vn2tc, ott_graph)
    diff_name_matched.extend(dnm)
    out.write('{} in {} and {}:\n'.format(len(both), ott_id, res_id))
    both.extend(diff_name_matched)
    _write_report_info_matched(out, both, None, None)


def _write_report_info_unmatched(out, just_in, obj_lookup, other_graph):
    n = 1
    diff_name_matched = []
    for name in just_in:
        found_tc = obj_lookup[name]
        to_extend = []
        if found_tc.is_specimen_based:
            gn = found_tc.has_name.genus_name
            sn = found_tc.has_name.sp_epithet
            if gn and sn:
                potential_genera = other_graph.find_valid_genus(gn.name)
                gen_with_correct_epi = []
                for genus in potential_genera:
                    gen_with_correct_epi.extend(genus.find_valid_species(gn.name, sn.name))
                to_extend.extend(gen_with_correct_epi)
        if to_extend:
            other_n = [i.canonical_name for i in to_extend]
            diff_name_matched.append((name, other_n, found_tc, to_extend))
        else:
            out.write('{} "{}" : '.format(n, name))
            if obj_lookup is not None:
                obj = obj_lookup[name]
                obj.explain(out)
            out.write('\n')
            n += 1
    return diff_name_matched

def _write_report_info_matched(out, just_in, obj_lookup, other_obj_lookup):
    for n, i in enumerate(just_in):
        alnl = i[1]
        if i[0] == alnl:
            o = ''
        else:
            o = ' => "{}"'.format('", "'.join([j.name for j in alnl]))
        out.write('{} "{}" {} '.format(1 + n, i[0], o))
        out.write('\n')


