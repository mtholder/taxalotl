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
            both.append((ott_name, tax_con, ott_name, ref_tc))
        else:
            just_ott.append((ott_name, tax_con))
    for ref_name, tax_con in ref_vn2tc.items():
        if ref_name not in ott_vn2tc:
            just_ref.append((ref_name, tax_con))
    just_ott.sort()
    just_ref.sort()
    dnm, ott_unmatched = _gather_unmatched(just_ott, ref_graph)
    diff_name_matched = dnm
    dnm, ref_unmatched = _gather_unmatched(just_ref, ott_graph, rev_order=True)
    diff_name_matched.extend(dnm)
    both.extend(diff_name_matched)
    ott_matched_set, ref_matched_set = set(), set()
    for el in both:
        ott_tc, ref_tc = el[1], el[3]
        ott_matched_set.add(ott_tc)
        if ref_tc is None:
            continue
        m = ref_matched_set.update if isinstance(ref_tc, list) else ref_matched_set.add
        m(ref_tc)
    ott_unmatched = [i for i in ott_unmatched if i[1] not in ott_matched_set]
    ref_unmatched = [i for i in ref_unmatched if i[1] not in ref_matched_set]
    # Report
    out.write('{} only in  {} :\n'.format(len(just_ott), ott_id))
    _write_tc_status(out, ott_unmatched)
    out.write('{} only in  {} :\n'.format(len(just_ref), res_id))
    _write_tc_status(out, ref_unmatched)
    out.write('{} in {} and {}:\n'.format(len(both), ott_id, res_id))
    _write_report_info_matched(out, both)
    out.write('Checking type of names of synonyms for {} \n'.format(ott_id))
    _check_syn_name_types(out, ott_graph)
    out.write('Checking type of names of synonyms for {} \n'.format(res_id))
    _check_syn_name_types(out, ref_graph)

def _check_syn_name_types(out, graph):
    msg_tmp = '"{}" is{} specimen-typed, but its synonym "{}" is{}.\n'
    for valid_tc in graph.taxon_concept_list:
        if not valid_tc.is_the_valid_name:
            continue
        oisb = valid_tc.is_specimen_based
        vtcs, oths = ('', ' not') if oisb else (' not', '')
        for syn in valid_tc.synonym_list:
            if syn.is_specimen_based != oisb:
                vtn, othn = valid_tc.canonical_name.name, syn.canonical_name.name
                out.write(msg_tmp.format(vtn, vtcs, othn, oths))
        for syn in valid_tc.problematic_synonym_list:
            vtn, othn = valid_tc.canonical_name.name, syn['name']
            oths = ' flagged as "{}"'.format(syn['problem'])
            out.write(msg_tmp.format(vtn, vtcs, othn, oths))

def _gather_unmatched(just_in, other_graph, rev_order=False):
    n = 1
    diff_name_matched = []
    still_unmatched = []
    for name_tc in just_in:
        name, found_tc = name_tc
        to_extend = []
        if found_tc.is_specimen_based and not found_tc.hybrid:
            gn = found_tc.has_name.genus_name
            sn = found_tc.has_name.sp_epithet
            if gn and sn:
                potential_genera = other_graph.find_valid_genus(gn.name)
                gen_with_correct_epi = []
                if found_tc.rank == 'species':
                    for genus in potential_genera:
                        tta = []
                        if found_tc.undescribed:
                            tta = genus.find_undescribed_species_for_name(gn.name, sn.name)
                        else:
                            tta = genus.find_valid_species(gn.name, sn.name)
                        gen_with_correct_epi.extend(tta)
                to_extend.extend(gen_with_correct_epi)
        if to_extend:
            other_n = [i.canonical_name for i in to_extend]
            if rev_order:
                tup = (other_n, to_extend, name, found_tc)
            else:
                tup = (name, found_tc, other_n, to_extend)
            diff_name_matched.append(tup)
        else:
            still_unmatched.append(name_tc)
    return diff_name_matched, still_unmatched

def _write_report_info_matched(out, just_in):
    for n, i in enumerate(just_in):
        ott_name, ott_tc, ref_name, ref_tc_or_list = i
        if ott_name == ref_name:
            o = ''
        else:
            rtcl = ref_tc_or_list
            if isinstance(rtcl, list) and len(rtcl) == 1:
                rtcl = rtcl[0]
            if isinstance(rtcl, list):
                o = ' => ["{}"]'.format('", "'.join([j.canonical_name.name for j in rtcl]))
            else:
                o = ' => "{}"'.format(rtcl.canonical_name.name)
        out.write('{} "{}" {} '.format(1 + n, i[0], o))
        out.write('\n')


def _write_tc_status(out, name_tc_list):
    n = 1
    for blob in name_tc_list:
        name = blob[0]
        obj = blob[1]
        out.write('{} '.format(n))
        obj.explain(out)
        out.write('\n')
        n += 1