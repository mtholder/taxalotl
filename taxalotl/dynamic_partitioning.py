from __future__ import print_function

import os
import re
import copy

from peyotl import get_logger, read_as_json, write_as_json, assure_dir_exists

from taxalotl.partitions import (get_relative_dir_for_partition,
                                 get_all_taxdir_and_misc_uncles,
                                 )
from taxalotl.tax_partition import (TAX_SLICE_CACHE,
                                    get_taxon_partition,
                                    PartitionedTaxDirBase)

_LOG = get_logger(__name__)
_norm_char_pat = re.compile(r'[-a-zA-Z0-9._]')


def _escape_odd_char(s):
    l = []
    for i in s:
        if _norm_char_pat.match(i):
            l.append(i)
        else:
            l.append('_')
    return ''.join(l)


def perform_dynamic_separation(ott_res, part_key, sep_fn, suppress_cache_flush=False):
    """Called where part_key is a PART_NAME element from the OTT 3 separation taxa."""
    fragment = get_relative_dir_for_partition(part_key)
    general_dynamic_separation(ott_res,
                               fragment,
                               sep_fn,
                               suppress_cache_flush=suppress_cache_flush)


def general_dynamic_separation(ott_res,
                               fragment,
                               sep_fn,
                               suppress_cache_flush=False):
    """
    :param suppress_cache_flush: 
    :param ott_res: resource wrapper for OTT
    :param fragment: path_frag relative to the top of the partitioned dir
    :param sep_fn:  name of the separations file ('__sep__.json')
    :return: None
    """
    sep_json_filepath = os.path.join(ott_res.partitioned_filepath, fragment, sep_fn)
    if not os.path.isfile(sep_json_filepath):
        _LOG.info('No separators found for {}'.format(fragment))
        return
    sep_obj = read_as_json(sep_json_filepath)
    try:
        _general_dynamic_separation_from_obj(ott_res,
                                             fragment,
                                             sep_obj,
                                             sep_fn)
    finally:
        if not suppress_cache_flush:
            TAX_SLICE_CACHE.flush()


class VirtualTaxonomyToRootSlice(PartitionedTaxDirBase):
    """Represents a taxon for a source, and all of "uncles" back to the root of the taxonomy.
    """

    def __init__(self,
                 res,
                 fragment):
        _LOG.info("Creating VirtualTaxonomyToRootSlice for {}".format(fragment))
        PartitionedTaxDirBase.__init__(self, res, fragment)
        self._taxon_partition = None
        self._garbage_bin_tax_part = None
        self._to_remove = None
        self.misc_uncle = None
        self._has_flushed = False
        # find misc_uncle
        ad = get_all_taxdir_and_misc_uncles(res.partitioned_filepath, fragment, res.id)
        infile_name_list = [os.path.join(i, res.taxon_filename) for i in ad]
        assert infile_name_list
        term_most = infile_name_list.pop(0)
        assert term_most == self.tax_fp_unpartitioned
        if infile_name_list:
            uncle_fp = infile_name_list[0]
            uk = (VirtualTaxonomyToRootSlice,self.src_id, uncle_fp)
            self.misc_uncle = TAX_SLICE_CACHE.get(uk)
            if self.misc_uncle is None:
                uf = os.path.split(self.fragment)[0]
                mu = get_virtual_tax_to_root_slice(res,
                                                   fragment=uf)
                self.misc_uncle = mu

    def get_vttrs_for_fragment(self, fragment):
        ck = (VirtualTaxonomyToRootSlice, self.src_id, fragment)
        return TAX_SLICE_CACHE.get(ck)

    @property
    def taxon_partition(self):
        if self._taxon_partition is None:
            tp = get_taxon_partition(self.res, self.fragment)
            self._taxon_partition = tp
        return self._taxon_partition

    def separate(self,
                 fragment_to_partition,
                 list_of_subdirname_and_roots,
                 dest_tax_part_obj=None):
        assert not self._has_flushed
        if dest_tax_part_obj is None and (not list_of_subdirname_and_roots):
            m = "No {} mapping to separate {}"
            _LOG.info(m.format(self.src_id, fragment_to_partition))
            return
        if fragment_to_partition == self.fragment:
            assert dest_tax_part_obj is None
            dest_tax_part_obj = self.taxon_partition
            dest_tax_part_obj.do_partition(list_of_subdirname_and_roots)
        else:
            own_tp = self.taxon_partition
            assert dest_tax_part_obj is not None
            assert list_of_subdirname_and_roots is None
            assert own_tp
            assert own_tp is not dest_tax_part_obj
            own_tp.move_from_misc_to_new_part(dest_tax_part_obj)
        if self.misc_uncle is not None:
            self.misc_uncle.separate(fragment_to_partition,
                                     list_of_subdirname_and_roots=None,
                                     dest_tax_part_obj=dest_tax_part_obj)

    def flush(self):
        if self._has_flushed:
            return
        _LOG.info("flushing VirtualTaxonomyToRootSlice for {}".format(self.fragment))
        if self._taxon_partition:
            tp = self._taxon_partition
            self._taxon_partition = None
            tp.flush()
        self._has_flushed = True
        TAX_SLICE_CACHE.try_del(self.cache_key)


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


def get_virtual_tax_to_root_slice(res,
                                  fragment):
    ck = (VirtualTaxonomyToRootSlice, res.id, fragment)
    c = TAX_SLICE_CACHE.get(ck)
    if c is not None:
        return c
    return VirtualTaxonomyToRootSlice(res, fragment)


def _general_dynamic_separation_from_obj(ott_res,
                                         fragment,
                                         sep_obj,
                                         sep_fn):
    """Separtes all sources and handles the recursion. Delegates to do_new_separation_of_src
    
    :param ott_res: 
    :param fragment: 
    :param sep_obj: 
    :param sep_fn: 
    :return: 
    """
    _LOG.info('breaking {} using {} separators: {}'.format(fragment,
                                                           len(sep_obj.keys()),
                                                           sep_obj.keys()))
    # Collect the set of source taxonomy IDs and mapping of sep ID (which is
    #   currently an OTT ID to the munged uniqname that will act as a readable
    #   subdirectory
    src_set = set()
    sep_id_to_fn = {}
    aug_sep_obj = {}
    for sep_id, i in sep_obj.items():
        a = copy.deepcopy(i)
        aug_sep_obj[sep_id] = a
        a['src_dict']["ott"] = [sep_id]
        src_set.update(a['src_dict'].keys())
        sep_id_to_fn[sep_id] = _escape_odd_char(i["uniqname"])
    _LOG.info('src_set: {}'.format(src_set))
    for src_id in src_set:
        if src_id != 'gbif':
            continue
        _gen_dyn_separation_from_obj_for_source(ott_res,
                                                fragment=fragment,
                                                sep_obj=sep_obj,
                                                src_id=src_id,
                                                sep_fn=sep_fn)


def _gen_dyn_separation_from_obj_for_source(ott_res,
                                            fragment,
                                            sep_obj,
                                            src_id,
                                            sep_fn):
    sep_id_to_fn = {}
    for sep_id, i in sep_obj.items():
        sep_id_to_fn[sep_id] = _escape_odd_char(i["uniqname"])
    taxalotl_config = ott_res.config
    res = taxalotl_config.get_terminalized_res_by_id(src_id)
    virt_taxon_slice = get_virtual_tax_to_root_slice(res, fragment)
    try:
        sep_dir_ids_list = []
        for sep_id, i in sep_obj.items():
            t = (sep_id_to_fn[sep_id], i['src_dict'][src_id])
            sep_dir_ids_list.append(t)
        _LOG.info("frag {} sep_dir_ids_list={}".format(fragment, sep_dir_ids_list))
        virt_taxon_slice.separate(fragment, sep_dir_ids_list)

        # Write the subdirectory separtion files, if needed (only first source)
        # and gather the list of recursive calls
        recursive_call_list = []
        for sep_id, obj in sep_obj.items():
            sub = obj.get("sub")
            if not sub:
                continue
            next_frag = os.path.join(fragment, sep_id_to_fn[sep_id])
            subdir = os.path.join(ott_res.partitioned_filepath, next_frag)
            assure_dir_exists(subdir)
            subjson = os.path.join(subdir, sep_fn)
            if not os.path.exists(subjson):
                write_as_json(sub, subjson, indent=2, sort_keys=True, separators=(',', ': '))
                _LOG.info('sub separation file written to "{}"'.format(subjson))
            rec_el = (next_frag, sub)
            recursive_call_list.append(rec_el)
        #_LOG.info("Skipping recursion"); recursive_call_list = []
        for recurse_el in recursive_call_list:
            next_frag, next_sep_obj = recurse_el
            _gen_dyn_separation_from_obj_for_source(ott_res=ott_res,
                                                    fragment=next_frag,
                                                    sep_obj=next_sep_obj,
                                                    src_id=src_id,
                                                    sep_fn=sep_fn)
    finally:
        virt_taxon_slice.flush()
