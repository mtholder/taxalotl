#!/usr/bin/env python
# from __future__ import print_function

import os

from peyotl import get_logger
from peyotl.utility.str_util import StringIO

from taxalotl.ott_schema import (read_taxonomy_to_get_single_taxon,
                                 )
from taxalotl.partitions import (fill_empty_anc_of_mapping,
                                 get_fragment_from_part_name,
                                 MISC_DIRNAME,
                                 PREORDER_PART_LIST,
                                 PARTS_BY_NAME)
from taxalotl.resource_wrapper import ResourceWrapper, TaxonomyWrapper
from taxalotl.tax_partition import (get_taxon_partition, get_root_ids_for_subset, )

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


def ott_fetch_root_taxon_for_partition(res, parts_key, root_id):
    tax_dir = res.get_taxdir_for_root_of_part(parts_key)
    if not tax_dir:
        _LOG.info('No taxon file found for {}'.format(parts_key))
        return None
    # _LOG.info('{} root should be in {}'.format(parts_key, tax_dir))
    taxon_obj = read_taxonomy_to_get_single_taxon(tax_dir, root_id)
    if not tax_dir:
        _LOG.info(
            'Root taxon for {} with ID {} not found in {}'.format(parts_key, root_id, tax_dir))
        return None
    return taxon_obj


def ott_build_paritition_maps(res):
    # srcs = {i[0]:i[1] for i in res.inputs if i[0] and is_str_type(i[1])}
    src_pre_to_map = {}
    all_part_keys = set()
    for part_key in PREORDER_PART_LIST:
        if part_key != 'Life':
            all_part_keys.add(part_key)
        all_part_keys.update(PARTS_BY_NAME[part_key])
    for pk in all_part_keys:
        if pk == MISC_DIRNAME:
            continue
        tax_dir = res.get_taxdir_for_part(pk)
        rids = get_root_ids_for_subset(tax_dir)
        _LOG.info("{} -> {} roots = {}".format(pk, tax_dir, rids))
        if not rids:
            continue
        assert len(rids) == 1
        root_id = rids.pop()
        taxon = ott_fetch_root_taxon_for_partition(res, pk, root_id)
        _LOG.info('root for {} has sources: {}'.format(pk, taxon.src_dict))
        for src_pre, src_ids in taxon.src_dict.items():
            src_pre_to_map.setdefault(src_pre, {})[pk] = src_ids
    import copy
    filled = {}
    for src_pre, mapping in src_pre_to_map.items():
        pm = copy.deepcopy(mapping)
        fill_empty_anc_of_mapping(mapping)
        filled[src_pre] = {k: list(v) for k, v in mapping.items()}
        for k, v in mapping.items():
            pv = pm.get(k)
            if pv is None:
                _LOG.info("{}: {} -> (None -> {})".format(src_pre, k, v))
            else:
                assert v == pv
    if 'silva' in filled:
        del filled['silva']
    # @TODO: make this automatic not generic
    gbif = filled['gbif']
    assert "Chloroplastida" not in gbif
    gbif["Chloroplastida"] = [13, 9, 7707728, 7819616, 36, 35]

    return filled


UNSTABLE_SRC_PREFIXES = frozenset(['h2007', 'study713', 'https', 'http'])


class NewSeparator(object):
    def __init__(self, ott_taxon_obj):
        self.taxon = ott_taxon_obj
        self.sub_separators = {}

    def __str__(self):
        out = StringIO()
        self.write_str(out)
        return out.getvalue()

    def write_str(self, out, pref=''):
        t = self.taxon
        if pref:
            pref = pref + '/'
        out.write('{}{} ott{} sources={}\n'.format(pref,
                                                   t.name_that_is_unique,
                                                   t.id,
                                                   t.src_dict))
        for name, sep in self.sub_separators.items():
            if pref:
                ns = pref + name
            else:
                ns = name
            sep.write_str(out, ns)

    def as_dict(self):
        d = {"name": self.taxon.name, "uniqname": self.taxon.name_that_is_unique}
        sd = {}
        for k, v in self.taxon.src_dict.items():
            vl = list(v)
            vl.sort()
            sd[k] = vl
        d["src_dict"] = sd
        s = {}
        for v in self.sub_separators.values():
            s[v.taxon.id] = v.as_dict()
        if s:
            d["sub"] = s
        return d

    def num_sub_separators(self):
        n = 0
        for el in self.sub_separators.values():
            n += 1 + el.num_sub_separators()
        return n


class NestedNewSeparator(object):
    def __init__(self):
        self.separators = {}

    def __bool__(self):
        return bool(self.separators)

    def add_separtors_for_tree(self, tree, sep_ids):
        self._add_separtors_for_descendants(self, tree.root, sep_ids, self.separators)

    def _add_separtors_for_descendants(self, nd, sep_ids, par_dict):
        child_refs = nd.children_refs
        for child in child_refs:
            if child.id in sep_ids:
                if nd.id in sep_ids and len(child_refs) == 1:
                    # suppress monotypic separators children as separators
                    self._add_separtors_for_descendants(child, sep_ids, par_dict)
                else:
                    ns = NewSeparator(child)
                    par_dict[child.id] = ns
                    self._add_separtors_for_descendants(child, sep_ids, ns.sub_separators)
            else:
                self._add_separtors_for_descendants(child, sep_ids, par_dict)

    def __str__(self):
        out = StringIO()
        for name, sep in self.separators.items():
            out.write('top-level name={}\n'.format(name))
            sep.write_str(out)
        return out.getvalue()

    def as_dict(self):
        d = {}
        for el in self.separators.values():
            d[el.taxon.id] = el.as_dict()
        return d

    def num_separators(self):
        n = 0
        for el in self.separators.values():
            n += 1 + el.num_sub_separators()
        return n


def _add_nst_subtree_el_to_dict(rd, nst_el, par_to_child):
    sep_taxon, children = nst_el
    if sep_taxon is not None:
        nst = NewSeparator(sep_taxon)
        nd = nst.sub_separators
        rd[sep_taxon.name_that_is_unique] = nst
    else:
        nd = rd
    for c in children:
        next_el = par_to_child[c]
        _add_nst_subtree_el_to_dict(nd, next_el, par_to_child)

NON_SEP_RANKS = frozenset(['forma', 'no rank - terminal', 'species',
                           'species group', 'species subgroup', 'varietas', 'variety', ])


def get_stable_source_keys(taxon):
    all_src_keys = taxon.src_dict.keys()
    filtered = [i for i in all_src_keys if not i.startswith('additions')]
    return [i for i in filtered if i not in UNSTABLE_SRC_PREFIXES]

# noinspection PyAbstractClass
class OTTaxonomyWrapper(TaxonomyWrapper):
    resource_type = 'open tree taxonomy'

    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    def diagnose_new_separators(self, current_partition_key):
        fragment = get_fragment_from_part_name(current_partition_key)
        tax_part = get_taxon_partition(self, fragment)
        tax_part.read_inputs_for_read_only()
        tax_forest = tax_part.get_taxa_as_forest()
        _LOG.info('{} taxon trees read from {}'.format(len(tax_forest.roots), tax_part.tax_fp))
        nns = NestedNewSeparator()
        for tree in tax_forest.trees:
            leaf_set = set()
            src_prefix_set = set()
            for v in tree.id_to_taxon.values():
                if v.children_refs is None:
                    leaf_set.add(v.id)
                filtered = get_stable_source_keys(v)
                src_prefix_set.update(filtered)
            max_num_srcs = len(src_prefix_set)
            _LOG.info("Relevant sources appear to be: {}".format(src_prefix_set))
            assert max_num_srcs > 0
            nst = set
            _LOG.info("root ids for current partition ({}) are: {}".format(current_partition_key, rids))
            for i, obj in tree.id_to_taxon.items():
                if obj is tree.root:
                    continue  # no point in partitioning at the root taxon
                if i in leaf_set:
                    continue  # no point in partitioning leaves...
                if len(get_stable_source_keys(obj)) == max_num_srcs:
                    if not obj.rank or (obj.rank not in NON_SEP_RANKS):
                        nst.add(i)
            nns.add_separtors_for_tree(tree, nst)
        if not nst:
            _LOG.debug('No new separators found for "{}"'.format(current_partition_key))
            return None
        rel_dir_for_part = get_fragment_from_part_name(current_partition_key)
        return {rel_dir_for_part: nns}

    def build_paritition_maps(self):
        return ott_build_paritition_maps(self)

    def get_primary_partition_map(self):
        return OTT_PARTMAP


# noinspection PyAbstractClass
class OTTaxonomyIdListWrapper(ResourceWrapper):
    resource_type = 'open tree taxonomy idlist'

    def __init__(self, obj, parent=None, refs=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs)

    def has_been_normalized(self):
        dfp = self.normalized_filepath
        return (dfp is not None
                and os.path.exists(dfp)
                and os.path.exists(os.path.join(dfp, 'by_qid.csv')))
