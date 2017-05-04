#!/usr/bin/env python
from __future__ import print_function
from taxalotl.interim_taxonomy_struct import InterimTaxonomyData
from peyotl.phylo.tree import parse_newick
from peyotl import get_logger
import os

_LOG = get_logger(__name__)


def normalize_newick(unpacked_fp, normalized_fp, resource_wrapper):
    lfp = os.path.join(unpacked_fp, resource_wrapper.local_filename)
    x = parse_newick(filepath=lfp)
    itd = InterimTaxonomyData()
    for taxon_id, nd in enumerate(x.preorder_node_iter()):
        nd.taxon_id = taxon_id
        par = nd.parent
        label = nd.id
        if not label:
            raise RuntimeError("Expecting every node to have a label in an input newick")
        if '=' in label:
            rank, name = [i.strip() for i in label.split('=')]
            itd.to_rank[taxon_id] = rank
        else:
            name = label
        itd.register_id_and_name(taxon_id, name)
        if par is None:
            itd.to_par[taxon_id] = None
            itd.root_nodes.add(taxon_id)
        else:
            itd.to_par[taxon_id] = par.taxon_id
            # print(taxon_id, name, nd.parent)
    for nd in x.preorder_node_iter():
        children = nd.children
        if children:
            itd.to_children[nd.taxon_id] = [c.taxon_id for c in children]
    itd.write_to_dir(normalized_fp)
