from __future__ import print_function

import codecs
from peyotl import get_logger
from taxalotl.partitions import do_partition, separate_part_list, TaxAndSynFileOnlyPartitionElement
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

