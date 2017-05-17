from __future__ import print_function

import codecs
from peyotl import get_logger
from taxalotl.partitions import do_partition, separate_part_list, TaxAndSynFileOnlyPartitionElement
_LOG = get_logger(__name__)

OTT_PARTMAP = {
    'Archaea': frozenset([996421]),
    'Bacteria': frozenset([844192]),
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

OTT_3_SEPARATION_TAXA = OTT_PARTMAP

def partition_ott(parts_dir, wrapper, part_name, part_keys, par_frag):
    do_partition(parts_dir,
                     wrapper,
                     part_name,
                     part_keys,
                     par_frag,
                     pe_class=TaxAndSynFileOnlyPartitionElement,
                     taxon_filename='taxonomy.tsv',
                     master_map=OTT_PARTMAP,
                     parse_and_partition_fn=_partition_ott_by_root_id
                     )

def _partition_col_by_root_id(complete_fp, partition_el_list):
    roots_set, by_roots, garbage_bin = separate_part_list(partition_el_list)
    id_to_line = {}
    id_by_par = {}
    syn_by_id = {}
    id_to_el = {}
    with codecs.open(complete_fp, 'rU', encoding='utf-8') as inp:
        iinp = iter(inp)
        header = iinp.next()
        prev_line = None
        # vt = unicode('\x0b') # Do some lines have vertical tabs? Of course they do....
        # istwo = unicode('\x1e')
        for n, line in enumerate(iinp):
            if not line.endswith('\n'):
                if prev_line:
                    prev_line = prev_line + line[:-1]
                else:
                    prev_line = line[:-1]
                continue
            elif prev_line:
                line = prev_line + line
                prev_line = ''
            ls = line.split('\t')
            if n % 1000 == 0:
                _LOG.info(' read line {}'.format(n))
            try:
                col_id, accept_id, par_id = ls[0], ls[4], ls[5]
                col_id = int(col_id)
                if accept_id:
                    accept_id = int(accept_id)
                    syn_by_id.setdefault(accept_id, []).append((col_id, line))
                else:
                    if col_id in roots_set:
                        match_l = [i[1] for i in by_roots if col_id in i[0]]
                        assert len(match_l) == 1
                        match_el = match_l[0]
                        id_to_el[col_id] = match_el
                        match_el.add(col_id, line)
                        if garbage_bin is not None:
                            garbage_bin.add(col_id, line)
                    else:
                        if par_id:
                            par_id = int(par_id)
                        match_el = id_to_el.get(par_id)
                        if match_el is not None:
                            id_to_el[col_id] = match_el
                            match_el.add(col_id, line)
                        else:
                            id_by_par.setdefault(par_id, []).append(col_id)
                            id_to_line[col_id] = line
            except:
                _LOG.exception("Exception parsing line {}:\n{}".format(n, line))
                raise
    while True:
        par_id_matched = [p for p in id_by_par.keys() if p in id_to_el]
        if not par_id_matched:
            break
        for p in par_id_matched:
            id_list = id_by_par[p]
            match_el = id_to_el[p]
            for col_id in id_list:
                id_to_el[col_id] = match_el
                match_el.add(col_id, id_to_line[col_id])
                del id_to_line[col_id]
            del id_by_par[p]
    if garbage_bin is not None:
        for col_id, line in id_to_line.items():
            garbage_bin.add(col_id, line)
    for accept_id, i_l_list in syn_by_id.items():
        match_el = id_to_el.get(accept_id)
        if match_el is None:
            match_el = garbage_bin
        for col_id, line in i_l_list:
            match_el.add_synonym(col_id, line)
    for part in partition_el_list:
        part.write_lines(header)
        pr = [r for r in roots_set if id_to_el.get(r) is part]
        part.write_roots(pr)

# Unused separation taxa:
# cellular organisms	93302
# Eukaryota	304358

