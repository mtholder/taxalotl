#!/usr/bin/env python
# from __future__ import print_function

import os
from collections import defaultdict

from peyotl import get_logger
from peyotl.utility.str_util import StringIO

from taxalotl.ott_schema import (read_taxonomy_to_get_single_taxon,
                                 )
from taxalotl.partitions import (fill_empty_anc_of_mapping,
                                 MISC_DIRNAME,
                                 PREORDER_PART_LIST,
                                 NAME_TO_PARTS_SUBSETS)
from taxalotl.resource_wrapper import ResourceWrapper, TaxonomyWrapper
from taxalotl.tax_partition import (get_roots_for_subset, )
from taxalotl.util import get_true_false_repsonse

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
        all_part_keys.update(NAME_TO_PARTS_SUBSETS[part_key])
    for pk in all_part_keys:
        if pk == MISC_DIRNAME:
            continue
        tax_dir = res.get_taxdir_for_part(pk)
        roots = get_roots_for_subset(tax_dir, res.get_misc_taxon_dir_for_part(pk))
        _LOG.info("{} -> {} roots = {}".format(pk, tax_dir, roots))
        if not roots:
            continue
        assert len(roots) == 1
        root_id, taxon = roots.popitem()
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
    # @TODO: hard copying the fact that OTT was built with SILVA-115
    if 'silva' in filled:
        filled['silva_115'] = filled['silva']
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
        self._add_separtors_for_descendants(tree.root, sep_ids, self.separators)

    def _add_separtors_for_descendants(self, nd, sep_ids, par_dict):
        child_refs = nd.children_refs
        if not child_refs:
            return
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
MIN_SEP_SIZE = 20


def get_stable_source_keys(taxon):
    all_src_keys = taxon.src_dict.keys()
    filtered = [i for i in all_src_keys if not i.startswith('additions')]
    r = [i for i in filtered if i not in UNSTABLE_SRC_PREFIXES]
    r.sort()
    return r


def _diagnose_relevant_sources(tree):
    source_to_count = defaultdict(int)
    leaf_set = set()
    src_prefix_set = set()
    for v in tree.id_to_taxon.values():
        if v.children_refs is None:
            leaf_set.add(v.id)
        filtered = get_stable_source_keys(v)
        src_prefix_set.update(filtered)
        for src in filtered:
            source_to_count[src] += 1
    # We'll define a "minor" source as one that is found in 5000 fold fewer taxa as
    #   the most common source for a slice. This is arbitrary.
    FOLD_DIFF = 5000
    max_seen = max(source_to_count.values())
    min_cutoff = 0 if max_seen < FOLD_DIFF else max_seen // FOLD_DIFF
    ac_src = [k for k, v in source_to_count.items() if v >= min_cutoff]
    ac_src.sort()
    _LOG.info("Relevant sources appear to be: {}".format(ac_src))
    return frozenset(ac_src)


DEFAULT_REL_SRC_SET = frozenset(['gbif', 'irmng', 'ncbi', 'worms'])
NO_WORMS_REL_SRC_SET = frozenset(['gbif', 'irmng', 'ncbi'])
SILVA_NOT_WORMS_REL_SRC_SET = frozenset(['gbif', 'irmng', 'ncbi', 'silva'])
PART_KEY_TO_REL_SRC_SET = { 'Insecta': NO_WORMS_REL_SRC_SET,
                            'Hymenoptera': NO_WORMS_REL_SRC_SET,
                            'Diptera': NO_WORMS_REL_SRC_SET,
                            'Coleoptera': NO_WORMS_REL_SRC_SET,
                            'Lepidoptera': NO_WORMS_REL_SRC_SET,
                            'Chordata': NO_WORMS_REL_SRC_SET,
                            'Fungi': NO_WORMS_REL_SRC_SET,
                            'Bacteria': SILVA_NOT_WORMS_REL_SRC_SET,
                            'Archaea': SILVA_NOT_WORMS_REL_SRC_SET,
                            }


def add_confirmed_sep(nns, tree, list_num_id_taxon, sep_name):
    if list_num_id_taxon:
        r = tree.root
        m = 'The current partition subtree "{}" has {} tips below it.'
        _LOG.info(m.format(r.name_that_is_unique, r.num_tips_below))
    top_sep_set = set()
    for nt, i, obj in list_num_id_taxon:
        if top_sep_set:
            already_sep = False
            for anc_obj in tree.to_root_gen(obj):
                if anc_obj.id in top_sep_set:
                    already_sep = True
                    break
            if already_sep:
                continue
        if sep_name is None:
            m = '"{}" has {} tips below it.'.format(obj.name_that_is_unique, nt)
            p = '{} Enter (y) to treat is a separator: '.format(m)
            if get_true_false_repsonse(p, def_value=True):
                top_sep_set.add(i)
        else:
            top_sep_set.add(i)
    if top_sep_set:
        nns.add_separtors_for_tree(tree, top_sep_set)


# noinspection PyAbstractClass
class OTTaxonomyWrapper(TaxonomyWrapper):
    resource_type = 'open tree taxonomy'
    schema = {'ott'}

    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    def get_source_for_sep_or_part(self, current_partition_key):
        try:
            return PART_KEY_TO_REL_SRC_SET.get(current_partition_key, DEFAULT_REL_SRC_SET)
        except:
            dn = self.config.get_fragment_from_part_name(current_partition_key)
            higher = os.path.split(os.path.split(dn)[0])[-1]
            if higher.lower() == 'life':
                raise ValueError('partition/separator name "{}" not found'.format(current_partition_key))
            return self.get_source_for_sep_or_part(higher)

    def diagnose_new_separators(self, current_partition_key, sep_name):
        tax_forest = self.get_taxon_forest_for_partition(current_partition_key)
        nns = NestedNewSeparator()
        if sep_name is None:
            lsn, lsep = 0, None
        else:
            lsn, lsep = len(sep_name), sep_name.lower()
        if tax_forest:
            try:
                ac_src = self.get_source_for_sep_or_part(current_partition_key)
            except:
                m = 'Might try calling _diagnose_relevant_sources for "{}", but in Jan 2019 moved to ' \
                    'making this hard coded in PART_KEY_TO_REL_SRC_SET.'
                raise NotImplementedError(m.format(current_partition_key))
            _LOG.info("Relevant sources for {} are recorded as to be: {}. size = {}".format(current_partition_key,
                                                                                 ac_src, len(ac_src)))
            for tree in tax_forest.trees:
                tree.add_num_tips_below()
                assert ac_src
                nst = []
                for i, obj in tree.id_to_taxon.items():
                    if obj is tree.root:
                        continue  # no point in partitioning at the root taxon
                    if not obj.children_refs:
                        continue  # no point in partitioning leaves...
                    if sep_name is not None:
                        obn, obun = obj.name, obj.uniqname
                        if (obn and (len(obn) == lsn and obn.lower() == lsep)) \
                            or (obun and (len(obun) == lsn and obun.lower() == lsep)):
                            nst.append((obj.num_tips_below, i, obj))
                            break
                    else:
                        sk_for_obj = set(get_stable_source_keys(obj))
                        if sk_for_obj.issuperset(ac_src):
                            if obj.num_tips_below >= MIN_SEP_SIZE:
                                if not obj.rank or (obj.rank not in NON_SEP_RANKS):
                                    nst.append((obj.num_tips_below, i, obj))
                nst.sort(reverse=True)
                add_confirmed_sep(nns, tree, nst, sep_name)
        if len(nns.separators) == 0:
            _LOG.info('No new separators found for "{}"'.format(current_partition_key))
            return None
        rel_dir_for_part = self.config.get_fragment_from_part_name(current_partition_key)
        return {rel_dir_for_part: nns}

    def build_paritition_maps(self):
        return ott_build_paritition_maps(self)

    def get_primary_partition_map(self):
        return OTT_PARTMAP

    def input_dict(self):
        d = {}
        for k, v in self.inputs:
            d[k] = v
        return d


# noinspection PyAbstractClass
class OTTaxonomyIdListWrapper(ResourceWrapper):
    resource_type = 'open tree taxonomy idlist'
    schema = {resource_type, 'ott id csv'}
    _norm_filename = 'by_qid.csv'

    def __init__(self, obj, parent=None, refs=None):
        ResourceWrapper.__init__(self, obj, parent=parent, refs=refs)
