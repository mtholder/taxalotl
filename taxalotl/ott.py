from __future__ import print_function

import codecs
from peyotl import get_logger
from taxalotl.partitions import (do_partition,
                                 separate_part_list,
                                 get_root_ids_for_subset,
                                 get_relative_dir_for_partition)
from taxalotl.interim_taxonomy_struct import read_taxonomy_to_get_id_to_fields

_LOG = get_logger(__name__)

OTT_PARTMAP = {
    'Archaea': frozenset([996421]),
    'Bacteria': frozenset([844192]),
    'Eukaryota': frozenset([304358]),
    'SAR': frozenset([5246039]),
    'Haptophyta': frozenset([151014]),
    'Rhodophyta': frozenset([878953]),
    'Archaeplastida': frozenset([5268475]),
    'Glaucophyta': frozenset([664970]),
    'Chloroplastida': frozenset([361838]),
    'Fungi': frozenset([352914]),
    'Metazoa': frozenset([691846]),
    'Annelida': frozenset([941620]),
    'Arthropoda': frozenset([632179]),
    'Malacostraca': frozenset([212701]),
    'Arachnida': frozenset([511967]),
    'Insecta': frozenset([1062253]),
    'Diptera': frozenset([661378]),
    'Coleoptera': frozenset([865243]),
    'Lepidoptera': frozenset([965954]),
    'Hymenoptera': frozenset([753726]),
    'Bryozoa': frozenset([442934]),
    'Chordata': frozenset([125642]),
    'Cnidaria': frozenset([641033]),
    'Ctenophora': frozenset([641212]),
    'Mollusca': frozenset([802117]),
    'Nematoda': frozenset([395057]),
    'Platyhelminthes': frozenset([555379]),
    'Porifera': frozenset([67819]),
    'Viruses': frozenset([4807313]),
}

# Unused separation taxa: cellular organisms	93302


OTT_3_SEPARATION_TAXA = OTT_PARTMAP


def partition_ott(res_wrapper, part_name, part_keys, par_frag):
    do_partition(res_wrapper,
                 part_name,
                 part_keys,
                 par_frag,
                 master_map=OTT_PARTMAP,
                 parse_and_partition_fn=_partition_ott_by_root_id)


def _partition_ott_by_root_id(complete_taxon_fp, syn_fp, partition_el_list):
    roots_set, by_roots, garbage_bin = separate_part_list(partition_el_list)
    id_to_line = {}
    id_by_par = {}
    syn_by_id = {}
    id_to_el = {}
    with codecs.open(syn_fp, 'rU', encoding='utf-8') as inp:
        iinp = iter(inp)
        syn_header = iinp.next()
        for n, line in enumerate(iinp):
            ls = line.split('\t|\t')
            if n % 1000 == 0:
                _LOG.info(' read synonym {}'.format(n))
            try:
                accept_id = int(ls[1])
                syn_by_id.setdefault(accept_id, []).append((None, line))
            except:
                _LOG.exception("Exception parsing line {}:\n{}".format(1 + n, line))
                raise
    with codecs.open(complete_taxon_fp, 'rU', encoding='utf-8') as inp:
        iinp = iter(inp)
        header = iinp.next()
        for n, line in enumerate(iinp):
            ls = line.split('\t|\t')
            if n % 1000 == 0:
                _LOG.info(' read taxon {}'.format(n))
            try:
                uid, par_id = ls[0], ls[1]
                uid = int(uid)
                if uid in roots_set:
                    match_l = [i[1] for i in by_roots if uid in i[0]]
                    assert len(match_l) == 1
                    match_el = match_l[0]
                    id_to_el[uid] = match_el
                    match_el.add(uid, line)
                    if garbage_bin is not None:
                        garbage_bin.add(uid, line)
                else:
                    if par_id:
                        par_id = int(par_id)
                    match_el = id_to_el.get(par_id)
                    if match_el is not None:
                        id_to_el[uid] = match_el
                        match_el.add(uid, line)
                    else:
                        id_by_par.setdefault(par_id, []).append(uid)
                        id_to_line[uid] = line
            except:
                _LOG.exception("Exception parsing line {}:\n{}".format(1 + n, line))
                raise
    return id_by_par, id_to_el, id_to_line, syn_by_id, roots_set, garbage_bin, header, syn_header


def ott_diagnose_new_separators(res, current_partition_key):
    tax_dir = res.get_taxdir_for_part(current_partition_key)
    rids = get_root_ids_for_subset(tax_dir)
    _LOG.info('tax_dir = {}'.format(tax_dir))
    id_to_obj = read_taxonomy_to_get_id_to_fields(tax_dir)
    _LOG.info('{} taxa read'.format(len(id_to_obj)))
    par_set = set()
    src_prefix_set = set()
    for v in id_to_obj.values():
        par_set.add(v.par_id)
        src_prefix_set.update(v.src_dict.keys())
    max_num_srcs = len(src_prefix_set)
    _LOG.info("Relevant sources appear to be: {}".format(src_prefix_set))
    nst = []
    if len(rids) > 1:
        rids = set()
    for i, obj in id_to_obj.items():
        if i in rids:
            continue  # no point in partitioning at the root taxon
        if i not in par_set:
            continue  # no point in partitioning leaves...
        if len(obj.src_dict) == max_num_srcs:
            nst.append((i, obj))
    if not nst:
        _LOG.debug('No new separators found for "{}"'.format(current_partition_key))
        return None
    par_to_child = {}
    to_par = {}
    for ott_id, obj in nst:
        par = obj.par_id
        to_par[ott_id] = par
        par_to_child.setdefault(par, [None, []])[1].append(ott_id)
        this_el = par_to_child.setdefault(ott_id, [None, []])
        assert this_el[0] is None
        this_el[0] = obj
    roots = set(par_to_child.keys()) - set(to_par.keys())
    rel_dir_for_part = get_relative_dir_for_partition(current_partition_key)
    return {rel_dir_for_part: NestedNewSeparator(roots, par_to_child)}


class NewSeparator(object):
    def __init__(self, ott_taxon_obj):
        self.taxon = ott_taxon_obj
        self.sub_separtors = {}


class NestedNewSeparator(object):
    def __init__(self, roots, par_to_child):
        ret_dict = {}
        for r in roots:
            curr_el = par_to_child[r]
            _add_nst_subtree_el_to_dict(ret_dict, curr_el, par_to_child)
        assert ret_dict
        self.separtors = ret_dict


def _add_nst_subtree_el_to_dict(rd, nst_el, par_to_child):
    sep_taxon, children = nst_el
    if sep_taxon is not None:
        nst = NewSeparator(sep_taxon)
        nd = nst.sub_separtors
        rd[sep_taxon.name_that_is_unique] = nst
    else:
        nd = rd
    for c in children:
        next_el = par_to_child[c]
        _add_nst_subtree_el_to_dict(nd, next_el, par_to_child)
