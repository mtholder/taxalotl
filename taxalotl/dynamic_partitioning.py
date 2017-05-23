from __future__ import print_function

import os
import re

from peyotl import get_logger, read_as_json, write_as_json

from taxalotl.partitions import (partition_from_mapping,
                                 get_relative_dir_for_partition,
                                 get_all_taxdir_and_misc_uncles,
                                 get_inverse_misc_non_misc_dir_for_tax)

_LOG = get_logger(__name__)

_norm_char_pat = re.compile(r'[-a-zA-Z0-9._]')


def escape_odd_char(s):
    l = []
    for i in s:
        if _norm_char_pat.match(i):
            l.append(i)
        else:
            l.append('_')
    return ''.join(l)


def perform_dynamic_separation(ott_res, part_key, sep_fn):
    """Called where part_key is a PART_NAME element from the OTT 3 separation taxa."""
    top_frag = get_relative_dir_for_partition(part_key)
    general_dynamic_separation(ott_res, top_frag, sep_fn)


class TaxonomySliceCache(object):
    def __init__(self):
        self._src_filepath_to_vttrs = {}

    def get(self, key):
        assert isinstance(key, tuple) and len(key) == 2
        return self._src_filepath_to_vttrs.get(key)

    def __setitem__(self, key, vttrs):
        assert isinstance(key, tuple) and len(key) == 2
        assert isinstance(vttrs, VirtualTaxonomyToRootSlice)
        old_val = self._src_filepath_to_vttrs.get(key)
        if old_val is not None and old_val is not vttrs:
            assert False, 'should not be creating a new VTTRS for an inp path!'
        self._src_filepath_to_vttrs[key] = vttrs

    def flush(self):
        kv = [(k, v) for k, v in self._src_filepath_to_vttrs.items()]
        _ex = None
        for k, v in kv:
            del self._src_filepath_to_vttrs[k]
            try:
                v.flush()
            except Exception as x:
                _ex = x
        if _ex is not None:
            raise _ex


def general_dynamic_separation(ott_res,
                               top_frag,
                               sep_fn,
                               tax_slice_cache=None):
    """
    :param ott_res: resource wrapper for OTT
    :param top_frag: path_frag relative to the top of the partitioned dir
    :param sep_fn:  name of the separations file ('__sep__.json')
    :param tax_slice_cache: 
    :return: None
    """
    sep_json_filepath = os.path.join(ott_res.partitioned_filepath, top_frag, sep_fn)
    if not os.path.isfile(sep_json_filepath):
        _LOG.info('No separators found for {}'.format(top_frag))
        return
    sep_obj = read_as_json(sep_json_filepath)
    cache_creator = False
    if tax_slice_cache is None:
        tax_slice_cache = TaxonomySliceCache()
        cache_creator = True
    try:
        _general_dynamic_separation_from_obj(ott_res,
                                             top_frag,
                                             sep_obj,
                                             sep_fn,
                                             tax_slice_cache=tax_slice_cache)
    except:
        if cache_creator:
            _LOG.exception("clearing TaxonomySliceCache")
            tax_slice_cache.flush()
        raise


class VirtualTaxonomyToRootSlice(object):
    def __init__(self, res, terminal_top_frag, tax_slice_cache, infile_name_list=None):
        self._taxon_partition = None
        self._to_remove = None
        self.res = res
        self.src_id = res.id
        pd = res.partitioned_filepath
        if infile_name_list is None:
            ad = get_all_taxdir_and_misc_uncles(pd, terminal_top_frag, res.id)
            infile_name_list = [os.path.join(i, res.taxon_filename) for i in ad]
            infile_name_list = [i for i in infile_name_list if os.path.isfile(i)]
            if not infile_name_list:
                raise ValueError("No input taxonomies for {} found in {}".format(res.id, ad))
        assert infile_name_list
        term_most = infile_name_list.pop(0)
        self.own_path = term_most
        self.inp_dir = os.path.split(term_most)[0]
        self.cache_key = (self.src_id, self.own_path)
        self.redundant_to = tax_slice_cache.get(self.cache_key)
        if self.redundant_to is not None:
            return
        if infile_name_list:
            uncle_fp = infile_name_list[0]
            self.misc_uncle = tax_slice_cache.get((self.src_id, uncle_fp))
            if self.misc_uncle is None:
                self.misc_uncle = get_virtual_tax_to_root_slice(res,
                                                                terminal_top_frag=None,
                                                                tax_slice_cache=tax_slice_cache,
                                                                infile_name_list=infile_name_list)
        else:
            self.misc_uncle = None
        self._has_flushed = False
        tax_slice_cache[self.cache_key] = self
        alt_dir, is_canon = get_inverse_misc_non_misc_dir_for_tax(self.inp_dir, res.id)
        alt_taxon_filepath = os.path.join(alt_dir, res.taxon_filepath)
        tax_slice_cache[(self.src_id, alt_taxon_filepath)] = self
        self.was_unpartitioned = is_canon
        self.tax_slice_cache = tax_slice_cache

    def separate(self, abs_par_dir, sep_dir_ids_list):
        assert self.redundant_to is None
        assert not self._has_flushed
        if not sep_dir_ids_list:
            _LOG.info("No {} mapping to separate {}".format(self.src_id, abs_par_dir))
            return
        if self._taxon_partition is None:
            tp, rm_file = partition_from_mapping(self.res,
                                                 sep_dir_ids_list,
                                                 inp_dir=self.inp_dir,
                                                 partition_parsing_fn=self.res.partition_parsing_fn,
                                                 par_dir=abs_par_dir)
            self._taxon_partition = tp
            assert self._to_remove is None
            self._to_remove = rm_file
        else:
            rm_file = self.repartition(sep_dir_ids_list, par_dir=abs_par_dir)
            if rm_file:
                assert self._to_remove is None or self._to_remove == rm_file
                self._to_remove = rm_file
        if self.misc_uncle is not None:
            self.misc_uncle.separate(abs_par_dir, sep_dir_ids_list)

    def flush(self):
        if self._taxon_partition:
            self._taxon_partition.write()
        if self._to_remove:
            _LOG.info("removing pre-partitioned file at {}".format(self._to_remove))
            os.unlink(self._to_remove)
        self._has_flushed = True

    def repartition(self, sep_dir_ids_list, par_dir):
        raise NotImplementedError("repartition")

"""
def do_new_separation_of_src(res,
                             abs_par_dir,
                             virt_taxon_slice,
                             sep_dir_and_id_sets_list):

    for inp_dir in src_inp_dirs:
        _LOG.info("partitioning new_par_dir = {} from {}".format(abs_par_dir, inp_dir))
        tp, rm_file = 
        tp_list.append(tp)
        to_remove.append(rm_file)
    merge_and_write_taxon_partition_list(tp_list)
    for to_remove_file in to_remove:
        if to_remove_file:
            _LOG.info("removing pre-partitioned file at {}".format(to_remove_file))
            os.unlink(to_remove_file)
"""


def get_virtual_tax_to_root_slice(res, terminal_top_frag, tax_slice_cache, infile_name_list=None):
    vttrs = VirtualTaxonomyToRootSlice(res,
                                       terminal_top_frag,
                                       tax_slice_cache,
                                       infile_name_list=infile_name_list)
    return vttrs if vttrs.redundant_to is None else vttrs.redundant_to


def _general_dynamic_separation_from_obj(ott_res,
                                         top_frag,
                                         sep_obj,
                                         sep_fn,
                                         tax_slice_cache):
    """Separtes all sources and handles the recursion. Delegates to do_new_separation_of_src
    
    :param ott_res: 
    :param top_frag: 
    :param sep_obj: 
    :param sep_fn: 
    :param tax_slice_cache: 
    :return: 
    """
    _LOG.info('breaking {} using {} separators: {}'.format(top_frag,
                                                           len(sep_obj.keys()),
                                                           sep_obj.keys()))
    # Collect the set of source taxonomy IDs and mapping of sep ID (which is
    #   currently an OTT ID to the munged uniqname that will act as a readable
    #   subdirectory
    src_set = set()
    sep_id_to_fn = {}
    for sep_id, i in sep_obj.items():
        src_set.update(i['src_dict'].keys())
        sep_id_to_fn[sep_id] = escape_odd_char(i["uniqname"])
    _LOG.info('src_set: {}'.format(src_set))
    #
    taxalotl_config = ott_res.config
    abs_par_dir = os.path.join(ott_res.partitioned_filepath, top_frag)
    for src_id in src_set:
        res = taxalotl_config.get_terminalized_res_by_id(src_id)
        virt_taxon_slice = VirtualTaxonomyToRootSlice(res, top_frag, tax_slice_cache)
        sep_dir_ids_list = []
        for sep_id, i in sep_obj.items():
            t = (sep_id_to_fn[sep_id], i['src_dict'][src_id])
            sep_dir_ids_list.append(t)
        virt_taxon_slice.separate(abs_par_dir, sep_dir_ids_list)
    # Write the subdirectory separtion files,
    recursive_call_list = []
    for sep_id, obj in sep_obj.items():
        sub = obj.get("sub")
        if not sub:
            continue
        next_frag = os.path.join(top_frag, sep_id_to_fn[sep_id])
        subdir = os.path.join(abs_par_dir, next_frag)
        subjson = os.path.join(subdir, sep_fn)
        write_as_json(sub, subjson, indent=2, sort_keys=True, separators=(',', ': '))
        _LOG.info('sub separation file written to "{}"'.format(subjson))
        rec_el = (next_frag, sub)
        recursive_call_list.append(rec_el)

    for recurse_el in recursive_call_list:
        next_frag, next_sep_obj = recurse_el
        _general_dynamic_separation_from_obj(ott_res=ott_res,
                                             top_frag=next_frag,
                                             sep_obj=next_sep_obj,
                                             sep_fn=sep_fn,
                                             tax_slice_cache=tax_slice_cache)
