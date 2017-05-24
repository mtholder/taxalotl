from __future__ import print_function

import os
import re

from peyotl import get_logger, read_as_json, write_as_json, assure_dir_exists

from taxalotl.partitions import (partition_from_mapping,
                                 get_relative_dir_for_partition,
                                 get_all_taxdir_and_misc_uncles,
                                 get_inverse_misc_non_misc_dir_for_tax,
                                 get_inp_taxdir,
                                 get_misc_inp_taxdir,)
from taxalotl.tax_partition import (TAX_SLICE_CACHE,
                                    get_taxon_partition,
                                    TaxonPartition,
                                    PartitionedTaxDirBase)

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

def general_dynamic_separation(ott_res,
                               top_frag,
                               sep_fn):
    """
    :param ott_res: resource wrapper for OTT
    :param top_frag: path_frag relative to the top of the partitioned dir
    :param sep_fn:  name of the separations file ('__sep__.json')
    :return: None
    """
    sep_json_filepath = os.path.join(ott_res.partitioned_filepath, top_frag, sep_fn)
    if not os.path.isfile(sep_json_filepath):
        _LOG.info('No separators found for {}'.format(top_frag))
        return
    sep_obj = read_as_json(sep_json_filepath)
    cache_creator = False
    try:
        _general_dynamic_separation_from_obj(ott_res,
                                             top_frag,
                                             sep_obj,
                                             sep_fn)
    except:
        TAX_SLICE_CACHE.flush()
        raise


class VirtualTaxonomyToRootSlice(PartitionedTaxDirBase):
    """Represents a taxon for a source, and all of "uncles" back to the root of the taxonomy.
    """
    def __init__(self,
                 res,
                 fragment):
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
            self.misc_uncle = TAX_SLICE_CACHE.get((self.src_id, uncle_fp))
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

    def separate(self, fragment_to_partition, sep_dir_ids_list):
        assert self.redundant_to is None
        assert not self._has_flushed
        if not sep_dir_ids_list:
            m = "No {} mapping to separate {}"
            _LOG.info(m.format(self.src_id, fragment_to_partition))
            return
        if fragment_to_partition == self.fragment:
            tp = self.taxon_partition
            tp.do_partition(sep_dir_ids_list)
            tp, rm_file = partition_from_mapping(self.res,
                                                 self.fragment,
                                                 sep_dir_ids_list,
                                                 inp_dir=self.inp_dir,
                                                 partition_parsing_fn=self.res.partition_parsing_fn,
                                                 par_dir=abs_par_dir)
            self._taxon_partition = tp
            assert self._to_remove is None
            self._to_remove = rm_file
            pd = self.res.partitioned_filepath + '/'
            for part in tp.partition_el:
                if part is tp.garbage_bin:
                    self._garbage_bin_tax_part = part
                else:
                    part_dir = part.get_taxon_partition_dir()
                    assert part_dir.startswith(pd)
                    part_frag = part_dir[len(pd):]
                    tf = part.taxon_fp
                    cached = TAX_SLICE_CACHE.get((self.src_id, tf))
                    if cached is None:
                        get_virtual_tax_to_root_slice(self.res,
                                                      terminal_top_frag=part_frag)
        else:
            rm_file = self.repartition(abs_par_dir, sep_dir_ids_list)
            if rm_file:
                assert self._to_remove is None or self._to_remove == rm_file
                self._to_remove = rm_file
        if self.misc_uncle is not None:
            self.misc_uncle.separate(fragment_to_partition, sep_dir_ids_list)

    def flush(self):
        if self._taxon_partition:
            self._taxon_partition.write()
        if self._to_remove and os.path.exists(self._to_remove):
            _LOG.info("removing pre-partitioned file at {}".format(self._to_remove))
            os.unlink(self._to_remove)
        self._has_flushed = True

    def repartition(self, par_dir, sep_dir_ids_list):
        assert self._garbage_bin_tax_part is not None
        _LOG.info('repartition: par_dir = {}'.format(par_dir))
        _LOG.info('repartition: sep_dir_ids_list = {}'.format(sep_dir_ids_list))
        _LOG.info('repartition: self.own_path = {}'.format(self.own_path))
        _LOG.info('repartition: garbage_bin.taxon_fp = {}'.format(self._garbage_bin_tax_part.taxon_fp))

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


def get_virtual_tax_to_root_slice(res,
                                  fragment):
    ck = (VirtualTaxonomyToRootSlice, res.id, fragment)
    c = TAX_SLICE_CACHE.get(ck)
    if c is not None:
        return c
    return VirtualTaxonomyToRootSlice(res, fragment)


def _general_dynamic_separation_from_obj(ott_res,
                                         top_frag,
                                         sep_obj,
                                         sep_fn):
    """Separtes all sources and handles the recursion. Delegates to do_new_separation_of_src
    
    :param ott_res: 
    :param top_frag: 
    :param sep_obj: 
    :param sep_fn: 
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
        if src_id != 'ncbi':
            continue
        res = taxalotl_config.get_terminalized_res_by_id(src_id)
        virt_taxon_slice = get_virtual_tax_to_root_slice(res, top_frag)
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
        subdir = os.path.join(ott_res.partitioned_filepath, next_frag)
        assure_dir_exists(subdir)
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
                                             sep_fn=sep_fn)
