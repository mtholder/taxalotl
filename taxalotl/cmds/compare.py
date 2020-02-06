#!/usr/bin/env python
from __future__ import print_function
from peyotl import get_logger
import sys
import os
from taxalotl.tax_partition import (get_taxon_partition, INP_TAXONOMY_DIRNAME, MISC_DIRNAME, OUTP_TAXONOMY_DIRNAME)
from taxalotl.util import OutFile, OutDir
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

def _init_bipartition_valid_taxa(ott_graph, ref_graph):
    ott_vn2tc = ott_graph.valid_name_to_taxon_concept_map
    ref_vn2tc = ref_graph.valid_name_to_taxon_concept_map
    just_ott, valid_in_both, just_ref = [], [], []
    for ott_name, tax_con in ott_vn2tc.items():
        if not tax_con.is_specimen_based:
            continue
        ref_tc = ref_vn2tc.get(ott_name)
        if ref_tc:
            valid_in_both.append((ott_name, tax_con, ott_name, [ref_tc]))
        else:
            just_ott.append((ott_name, tax_con))
    for ref_name, tax_con in ref_vn2tc.items():
        if ref_name not in ott_vn2tc:
            just_ref.append((ref_name, tax_con))
    just_ott.sort()
    just_ref.sort()
    ott_only_claims_valid, ott_unmatched = _gather_unmatched(just_ott, ref_graph)
    ref_only_claims_valid, ref_unmatched = _gather_unmatched(just_ref, ott_graph, rev_order=True)
    ott_matched_set, ref_matched_set = set(), set()
    for el in valid_in_both:
        ott_tc, ref_tc = el[1], el[3]
        ott_matched_set.add(ott_tc)
        if ref_tc is None:
            continue
        assert len(ref_tc) == 1
        ref_matched_set.update(ref_tc)
    ott_unmatched = [i for i in ott_unmatched if i[1] not in ott_matched_set]
    ref_unmatched = [i for i in ref_unmatched if i[1] not in ref_matched_set]
    bijection_with_1_invalid = []
    one_ott_to_many = []
    one_ref_to_many = []
    for el in ott_only_claims_valid:
        other_tc_list = el[-1]
        if len(other_tc_list) > 1:
            one_ott_to_many.append(el)
        else:
            bijection_with_1_invalid.append(el)
    if ref_only_claims_valid:
        other_tc_list = el[-1]
        if len(other_tc_list) > 1:
            one_ref_to_many.append(el)
        else:
            bijection_with_1_invalid.append(el)
    return {'valid_in_both': valid_in_both, 'bijection_with_1_invalid': bijection_with_1_invalid,
            'in_first_but_ambig': one_ott_to_many, 'in_second_but_ambig': one_ref_to_many,
            'found_only_in_first': ott_unmatched, 'valid_only_in_first': ott_only_claims_valid,
            'found_only_in_second': ref_unmatched, 'valid_only_in_second': ref_only_claims_valid,
            }

def _write_report(out, ott_id, ott_graph, res_id, ref_graph):
    # ott_vstc = ott_graph.valid_specimen_based_taxa
    # ref_vstc = ref_graph.valid_specimen_based_taxa
    out.write('comparing {} to {}\n'.format(ott_id, res_id))
    init_bipar = _init_bipartition_valid_taxa(ott_graph, ref_graph)
    # Report
    for key, id_to_show in [['found_only_in_first', ott_id], ['found_only_in_second', res_id]]:
        blob = init_bipar[key]
        out.write('{} only in  {} :\n'.format(len(blob), id_to_show))
        _write_tc_status(out, blob)
    for key, id_to_show in [['in_first_but_ambig', ott_id], ['in_second_but_ambig', res_id]]:
        blob = init_bipar[key]
        oid = ott_id if id_to_show == res_id else res_id
        out.write('{} in  {} mapping to >1 in {} :\n'.format(len(blob), id_to_show, oid))
        _write_tc_status(out, blob)
    out.write('Checking type of names of synonyms for {} \n'.format(ott_id))
    _check_syn_name_types(out, ott_graph)
    out.write('Checking type of names of synonyms for {} \n'.format(res_id))
    _check_syn_name_types(out, ref_graph)
    in_both = init_bipar['valid_in_both'] + init_bipar['bijection_with_1_invalid']
    _diagnose_higher_taxa(out, in_both, ott_id, ott_graph, res_id, ref_graph)

def _use_uniq_mapping_to_diagnose_higher_taxa(out, uniq_map, ott_id, ott_graph, res_id, ref_graph):
    ott2ref, ref2ott = _gen_uniq_tc_maps(uniq_map)
    _add_des_uniq_ids(ott_graph, ott2ref)
    _add_des_uniq_ids(ref_graph, ref2ott)
    ott_name_to_tc = ott_graph.canonical_name_str_to_taxon_concept_map
    ref_name_to_tc = ref_graph.canonical_name_str_to_taxon_concept_map
    checked_names = set()
    for ott_name, tax_con in ott_name_to_tc.items():
        if tax_con.is_specimen_based:
            continue
        checked_names.add(ott_name)
        if not tax_con.is_the_valid_name:
            rtc = ref_name_to_tc.get(ott_name)
            if rtc and rtc.is_the_valid_name:
                out.write('"{}" is valid in {}, but not in {}\n'.format(ott_name, res_id, ott_id))
            continue
        rtc = ref_name_to_tc.get(ott_name)
        if not rtc:
            out.write('"{}" is valid in {}, but not found in "{}"\n'.format(ott_name, ott_id, res_id))
            continue
        if not rtc.is_the_valid_name:
            out.write('"{}" is valid in {}, but not in "{}"\n'.format(ott_name, ott_id, res_id))
            continue
        proj_to_ref = set([ott2ref[i] for i in tax_con.des_uniq_ids])
        if proj_to_ref == rtc.des_uniq_ids:
            out.write('{} and {} agree on the definition of "{}"\n'.format(ott_id, res_id, ott_name))
        else:
            in_both = proj_to_ref.intersection(rtc.des_uniq_ids)
            missing_in_ott = rtc.des_uniq_ids - proj_to_ref
            missing_in_ref = proj_to_ref - rtc.des_uniq_ids
            missing_in_ott_str = '", "'.join([i.valid_name.name for i in missing_in_ott])
            missing_in_ref_str = '", "'.join([i.valid_name.name for i in missing_in_ref])
            mirs = ' {} contains ["{}"]'.format(ott_id, missing_in_ref_str) if missing_in_ref_str else ''
            mios = ' {} contains ["{}"].'.format(res_id, missing_in_ott_str) if missing_in_ott_str else '.'
            m = '"{}" shares {} uniq-mapping descendants. However:{}{}\n'
            out.write(m.format(ott_name, len(in_both), mirs, mios))

def _gen_uniq_tc_maps(uniq_map):
    ott2ref = {}
    ref2ott = {}
    for ott_tax_con in uniq_map:
        if isinstance(ott_tax_con, frozenset):
            assert len(ott_tax_con) == 1
            ott_tax_con = list(ott_tax_con)[0]
        assert isinstance(ott_tax_con.mapped_to, set)
        assert len(ott_tax_con.mapped_to) == 1
        ref_tax_con = next(iter(ott_tax_con.mapped_to))
        assert ott_tax_con not in ott2ref
        ott2ref[ott_tax_con] = ref_tax_con
        assert ref_tax_con not in ref2ott
        ref2ott[ref_tax_con] = ott_tax_con
        # sys.stdout.write('ott "{}" <-> cof "{}"\n'.format(ott_tax_con.valid_name.name, ref_tax_con.valid_name.name))
    return ott2ref, ref2ott

def _add_des_uniq_ids(graph, uniq_id_dict):
    for tax_con in graph.postorder_taxon_concepts():
        cs = tax_con.child_set
        tax_con.des_uniq_ids = set()
        if tax_con in uniq_id_dict:
            tax_con.des_uniq_ids.add(tax_con)
        if cs:
            for c in cs:
                tax_con.des_uniq_ids.update(c.des_uniq_ids)


def _set_mapped_to_and_partition(matched):
    # Reset the `mapped_to` field
    ott_tc_set, ref_tc_set = set(), set()
    for ott_name, ott_tc, ref_name, rtcl in matched:
        ott_tc.mapped_to = set()
        ott_tc_set.add(ott_tc)
        for rtc in rtcl:
            ref_tc_set.add(rtc)
            rtc.mapped_to = set()
    # initialize the `mapped_to` field
    for ott_name, ott_tc, ref_name, rtcl in matched:
        for rtc in rtcl:
            ott_tc.mapped_to.add(rtc)
            rtc.mapped_to.add(ott_tc)
    # partition into uniq_map, merge, split, snarl
    uniq_map, merge, split, snarl = set(), set(), set(), set()
    for ott_tc in ott_tc_set:
        np = len(ott_tc.mapped_to)
        if np == 1:
            ref_tc = next(iter(ott_tc.mapped_to))
            if len(ref_tc.mapped_to) == 1:
                uniq_map.add(ott_tc)
            else:
                assert len(ref_tc.mapped_to) > 1
                union = set()
                for i_ott_tc in ref_tc.mapped_to:
                    union.add(i_ott_tc)
                assert len(union) > 1
                key = frozenset(union)
                ref_union = set()
                dest = merge
                for i_ott_tc in key:
                    ref_union.update(i_ott_tc.mapped_to)
                    if len(ref_union) > 1:
                        dest = snarl
                        break
                dest.add(key)
        else:
            assert np > 1
            union = set()
            for ref_tc in ott_tc.mapped_to:
                union.update(ref_tc.mapped_to)
            is_snarl = len(union) > 1
            dest = snarl if is_snarl else split
            dest.add(frozenset(union))
    red_snarl = set()
    for x in snarl:
        add, remove = True, []
        for y in red_snarl:
            assert x != y
            if y.issubset(x):
                remove.append(y)
            if x.issubset(y):
                add = False
        if add:
            red_snarl.add(x)
        for td in remove:
            red_snarl.remove(td)
    return uniq_map, merge, split, red_snarl

def _diagnose_higher_taxa(out, matched, ott_id, ott_graph, res_id, ref_graph):
    uniq_map, merge, split, snarl = _set_mapped_to_and_partition(matched)
    # Report partitions
    out.write('{} names in {} need to be split into multiple names in and {}:\n'.format(len(split), ott_id, res_id))
    _write_using_mapped_to(out, split)
    out.write('{} sets of names in {} need to be merged into names in and {}:\n'.format(len(merge), ott_id, res_id))
    _write_using_mapped_to(out, merge)
    out.write('{} sets of names in {} need to be merged into names in a snarl in {}:\n'.format(len(snarl), ott_id, res_id))
    _write_using_mapped_to(out, snarl)
    out.write('{} names in {} map uniquely to a name in {}:\n'.format(len(uniq_map), ott_id, res_id))
    _write_using_mapped_to(out, uniq_map)
    # Use uniq_mapped to generate taxon_concept definitions
    _use_uniq_mapping_to_diagnose_higher_taxa(out, uniq_map, ott_id, ott_graph, res_id, ref_graph)


def _write_using_mapped_to(out, mapping_set):
    for n, k in enumerate(mapping_set):
        if isinstance(k, frozenset):
            union = set()
            for m in k:
                union.update(m.mapped_to)
            sl = [i.canonical_name.name for i in k]
            sl.sort()
            osl = '", "'.join(sl)
            usl = [i.canonical_name.name for i in union]
            usl.sort()
            if len(usl) == 1:
                out.write('{} ["{}"] -> "{}"\n'.format(1 + n, osl, usl[0]))
            else:
                out.write('{} ["{}"] -> ["{}"]\n'.format(1 + n, osl, '", "'.join(usl)))
        else:
            usl = [i.canonical_name.name for i in k.mapped_to]
            usl.sort()
            osl = k.canonical_name.name
            if len(usl) == 1:
                if osl == usl[0]:
                    out.write('{} "{}"\n'.format(1 + n, osl))
                else:
                    out.write('{} "{}" -> "{}"\n'.format(1 + n, osl, usl[0]))
            else:
                out.write('{} "{}" -> ["{}"]\n'.format(1 + n, osl, '", "'.join(usl)))


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


def _write_tc_status(out, name_tc_list):
    n = 1
    for blob in name_tc_list:
        if len(blob) == 2:
            obj = blob[1]
            out.write('  {} '.format(n))
            obj.explain(out)
            out.write('\n')
        else:
            matched_tc_list = blob[3]
            vnlist = [i.valid_name for i in matched_tc_list]
            o = '", "'.join([i.name for i in vnlist])
            out.write('"{}" ambiguous between ["{}"]\n'.format(blob[0], o))
        n += 1