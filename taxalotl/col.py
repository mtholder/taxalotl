from __future__ import print_function

import codecs
import os

from peyotl import (get_logger)

from taxalotl.partitions import (do_partition,
                                 TaxonFileOnlyPartitionElement,
                                 separate_part_list)

_LOG = get_logger(__name__)

COL_PARTMAP = {'Archaea': frozenset([33524792]),
               'Bacteria': frozenset([33521420]),
               'Eukaryota/__other__': frozenset([33521595, 33523363]),
               'Eukaryota/Archaeplastida': frozenset([33521293]),
               'Eukaryota/Fungi': frozenset([33521351]),
               'Eukaryota/Metazoa': frozenset([33521288]),
               'Eukaryota/Metazoa/Annelida': frozenset([33521477]),
               'Eukaryota/Metazoa/Arthropoda': frozenset([33521342]),
               'Eukaryota/Metazoa/Arthropoda/Malacostraca': frozenset([33521356]),
               'Eukaryota/Metazoa/Arthropoda/Arachnida': frozenset([33521366]),
               'Eukaryota/Metazoa/Arthropoda/Insecta': frozenset([33521474]),
               'Eukaryota/Metazoa/Arthropoda/Insecta/Diptera': frozenset([33521519]),
               'Eukaryota/Metazoa/Arthropoda/Insecta/Coleoptera': frozenset([33521475]),
               'Eukaryota/Metazoa/Arthropoda/Insecta/Lepidoptera': frozenset([33521695]),
               'Eukaryota/Metazoa/Arthropoda/Insecta/Hymenoptera': frozenset([33521644]),
               'Eukaryota/Metazoa/Bryozoa': frozenset([33524015]),
               'Eukaryota/Metazoa/Chordata': frozenset([33521289]),
               'Eukaryota/Metazoa/Cnidaria': frozenset([33522061]),
               'Eukaryota/Metazoa/Ctenophora': frozenset([33521313]),
               'Eukaryota/Metazoa/Mollusca': frozenset([33521301]),
               'Eukaryota/Metazoa/Nematoda': frozenset([33526516]),
               'Eukaryota/Metazoa/Platyhelminthes': frozenset([33521309]),
               'Eukaryota/Metazoa/Porifera': frozenset([33527549]),
               'Viruses': frozenset([33521407]),
               }


def partition_col(parts_dir, wrapper, part_name, part_keys, par_frag):
    do_partition(parts_dir,
                 wrapper,
                 part_name,
                 part_keys,
                 par_frag,
                 pe_class=TaxonFileOnlyPartitionElement,
                 taxon_filename='taxa.txt',
                 master_map=COL_PARTMAP,
                 parse_and_partition_fn=_partition_col_by_root_id
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
    return id_by_par, id_to_el, id_to_line, syn_by_id, roots_set, garbage_bin, header, None
