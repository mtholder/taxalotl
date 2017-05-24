#!/usr/bin/env python
# from __future__ import print_function
import codecs
import os

from peyotl import get_logger, assure_dir_exists

INP_TAXONOMY_DIRNAME = '__inputs__'
MISC_DIRNAME = '__misc__'
GEN_MAPPING_FILENAME = '__mapping__.json'

_LOG = get_logger(__name__)


class TaxonomySliceCache(object):
    def __init__(self):
        self._ck_to_obj = {}

    def get(self, key):
        assert isinstance(key, tuple) and len(key) == 3
        return self._ck_to_obj.get(key)


    def __setitem__(self, key, vttrs):
        assert isinstance(key, tuple) and len(key) == 3
        old_val = self._ck_to_obj.get(key)
        if old_val is not None and old_val is not vttrs:
            assert False, 'should not be creating a new object for a cached taxdi!'
        self._ck_to_obj[key] = vttrs

    def __getitem__(self, ck):
        return self._ck_to_obj[ck]

    def __delitem__(self, ck):
        obj = self._ck_to_obj.get(ck)
        if obj:
            del self._ck_to_obj[ck]
            obj.flush()

    def flush(self):
        kv = [(k, v) for k, v in self._ck_to_obj.items()]
        _ex = None
        for k, v in kv:
            del self._ck_to_obj[k]
            try:
                v.flush()
            except Exception as x:
                _LOG.exception('exception in flushing')
                _ex = x
        if _ex is not None:
            raise _ex
TAX_SLICE_CACHE = TaxonomySliceCache()

class PartitionedTaxDirBase(object):
    def __init__(self, res, fragment):
        self.fragment = fragment
        self.res = res
        self.src_id = res.id
        self.scaffold_dir = os.path.join(res.partitioned_filepath, fragment)
        self.tax_fp_unpartitioned = res.get_taxon_filepath_for_part(fragment)
        self.tax_fp_misc = res.get_misc_taxon_filepath_for_part(fragment)
        self.tax_dir_unpartitioned = res.get_taxon_dir_for_part(fragment)
        self.tax_dir_misc = res.get_misc_taxon_dir_for_part(fragment)
        sf = self.res.synonyms_filename
        self.syn_fp = os.path.join(self.tax_dir_unpartitioned, sf) if sf else None
        self.syn_fp_misc = os.path.join(self.tax_dir_misc, sf) if sf else None
        self.cache_key = (self.__class__, self.src_id, self.fragment)
        assert TAX_SLICE_CACHE.get(self.cache_key) is None
        TAX_SLICE_CACHE[self.cache_key] = self


    def scaffold_tax_subdirs(self):
        """Returns a list of subdirectory names for self.scaffold_dir with __misc__ suppressed"""
        n = []
        for x in os.listdir(self.scaffold_dir):
            if x == INP_TAXONOMY_DIRNAME or x == MISC_DIRNAME:
                continue
            if os.path.isdir(os.path.join(self.scaffold_dir, x)):
                n.append(x)
        return n

class LightTaxonomyHolder(object):
    def __init__(self, fragment):
        self.fragment = fragment
        self._id_order = []
        self._id_to_line = {}  # id -> line
        self._id_to_child_set = {}  # id -> set of child IDs
        self._syn_by_id = {}  # accepted_id -> list of synonym lines
        self._id_to_el = {}
        self._roots = set()
        self.taxon_header = None
        self.syn_header = None
        self.treat_syn_as_taxa = False

    def add_root_taxon_from_higher_tax_part(self, uid, par_id, line):
        self._roots(uid)
        self.add_taxon(uid, par_id, line)

    def add_taxon(self, uid, par_id, line):
        assert uid not in self._id_to_line
        self._id_to_line[uid] = line
        self._id_order.append(uid)
        self._id_to_child_set.setdefault(par_id, set()).add(uid)

    add_moved_taxon = add_taxon
    add_taxon_from_higher_tax_part = add_taxon

class PartitioningLightTaxHolder(LightTaxonomyHolder):
    def __init__(self, fragment):
        LightTaxonomyHolder.__init__(self, fragment)
        self._subdirname_to_tp_roots = {}
        self._misc_part = LightTaxonomyHolder()
        self._roots_for_sub = set()
        self._root_to_lth = {}

    def add_synonym(self, accept_id, syn_id, line):
        if self.treat_syn_as_taxa:
            # CoL uses the taxonomy file for synonyms.
            assert syn_id is not None
            self.read_taxon_line(syn_id, None, line)
        else:
            self._syn_by_id.setdefault(accept_id, []).append((syn_id, line))

    def read_taxon_line(self, uid, par_id, line):
        if uid in self._roots_for_sub:
            # _LOG.info("{} not in {}".format(uid, roots_set))
            match_el = self._root_to_lth[uid]
            self._id_to_el[uid] = match_el
            match_el.add_root_taxon_from_higher_tax_part(uid, par_id, line)
            # add as a leaf to the __misc__
            self._misc_part.add_moved_taxon(uid, par_id, line)
        else:
            if par_id:
                try:
                    par_id = int(par_id)
                except:
                    pass
            match_el = self._id_to_el.get(par_id)
            if match_el is not None:
                self._id_to_el[uid] = match_el
                match_el.add_taxon_from_higher_tax_part(uid, par_id, line)
            else:
                self._id_to_child_set.setdefault(par_id, set()).add(uid)
                self._id_to_line[uid] = line

    def _finish_partition_after_parse(self):
        """On entry _id_to_el will be set for the root elements (and some of their
            children), but taxa processed before their ancestors may have been missed.
        _id_to_child_list and _id_to_line are only filled for these
        """
        par_id_matched = [p for p in self._id_to_child_set.keys() if p in self._id_to_el]
        if not par_id_matched:
            break
        for p in par_id_matched:
            match_el = self._id_to_el[p]
            self._transfer_subtree(p, match_el)
        misc_el = self._misc_part
        mitcs = self._misc_part._id_to_child_set
        mitl = self._misc_part._id_to_line
        mite = self._misc_part._id_to_el
        for par_id, child_set in self._id_to_child_set.items():
            mics = mitcs.setdefault(par_id, set())
            mics.update(child_set)
            for c in child_set:
                curr_el = self._id_to_el.get(c)
                if curr_el:
                    assert False
                    # mite[c] = curr_el
                else:
                    self._id_to_el[c] = self._misc_part
                    line = self._id_to_line.get(c)
                    if not line:
                        m = "Child ID {} in {} not associated with a line or el"
                        raise ValueError(m.format(c, self.fragment))
                    mitl[c] = line
                    del self._id_to_line[c]

        for par_id in self._id_to_child_set.keys():
            if par_id not in mitl and par_id not in mite:
                assert par_id not in self._id_to_el
                line = self._id_to_line[par_id]
                if line:
                    self._id_to_el[par_id] = self._misc_part
                    mitl[par_id] = line
                    del self._id_to_line[par_id]
                else:
                    # par_id was mentioned by its children, but not included.
                    #   so each child is a root wrt its slice of the taxonomy
                    cset = self._id_to_child_set[par_id]
                    for c in cset:
                        dest_for_c = self._id_to_el[c]
                        dest_for_c._roots.update(cset)
        for k in self._id_to_line.keys():
            "End of partitioning and taxon {} has not been moved in {}"
            raise ValueError(m.format(k, self.fragment))
        # Save some memory
        self._id_to_child_set = None
        self._id_to_line = None

    def _transfer_subtree(self, par_id, part_element):
        child_list = self._id_to_child_set[par_id]
        self._id_to_el[par_id] = part_element
        line = self._id_to_line.get(par_id)
        if line is not None:
            part_element._id_to_line[par_id] = line
            del self._id_to_line[par_id]
        del self._id_to_child_set[par_id]
        for child_id in child_list:
            self._id_to_el[child_id] = part_element
            if child_id in self._id_to_child_set:
                self._transfer_subtree(child_id, part_element)
            else:
                part_element.add_taxon(child_id, par_id, self._id_to_line[child_id])
                del self._id_to_line[child_id]


    def _partition_synonyms(self):
        for accept_id, i_l_list in self.syn_by_id.items():
            match_el = self._id_to_el.get(accept_id)
            if match_el is None:
                m = "Transferring synonyms for {} to __misc__ of {}"
                _LOG.info(m.format(accept_id, self.fragment))
                match_el = self._misc_part
            for syn_id, line in i_l_list:
                match_el.add_synonym(accept_id, syn_id, line)
        # Save some memory
        self.syn_by_id = None
        self._id_to_el = None


class TaxonPartition(PartitionedTaxDirBase, PartitioningLightTaxHolder):
    def __init__(self, res, fragment):
        PartitioningLightTaxHolder.__init__(self, fragment)
        PartitionedTaxDirBase.__init__(self, res, fragment)
        self.treat_syn_as_taxa = self.syn_fp is None
        self._read = False
        self._read_from_misc = False
        self._fs_is_partitioned = None
        self._populated = False


    def _diagnose_state_of_fs(self):
        if os.path.exists(self.tax_fp_misc):
            self._fs_is_partitioned = True
        elif os.path.exists(self.tax_fp_unpartitioned):
            self._fs_is_partitioned = False
        else:
            self._fs_is_unpartitioned = None

    def do_partition(self, list_of_subdirname_and_roots):
        if self._subdirname_to_tp_roots:
            raise ValueError("do_partition called twice for {}".format(self.fragment))
        if not self._populated:
            self._diagnose_state_of_fs()
            if self._fs_is_partitioned is None:
                m = "Taxa files not found for {} and TaxonPartition is empty"
                raise ValueError(m.format(self.fragment))
        cur_sub_names = self.scaffold_tax_subdirs()
        if cur_sub_names:
            req_fulfilled = True
            for x in list_of_subdirname_and_roots:
                if x[0] not in cur_sub_names:
                    req_fulfilled = False
                    break
            if req_fulfilled:
                m = "Partition subdirs found for {} TaxonPartition"
                raise ValueError(m.format(self.fragment))
            m = "Some taxonomic subdirs found for {} TaxonPartition"
            raise ValueError(m.format(self.fragment))

        for subname, subroot in list_of_subdirname_and_roots:
            subfrag = os.path.join(self.fragment, subname)
            subtp = get_taxon_partition(self.res, subfrag)
            self._roots_for_sub.update(subroot)
            for r in subroot:
                self._root_to_lth[r] = subtp
            self._subdirname_to_tp_roots[subname] = (subtp, subroot)
        if self._populated:
            raise NotImplementedError("partition of populated via in-memory")
        else:
            self._read_inputs()

    def _read_inputs(self):
        self._read = True
        if os.path.exists(self.tax_fp_misc):
            self.tax_fp = self.tax_fp_misc
            self._read_from_misc = True
        else:
            self.tax_fp = self.tax_fp_unpartitioned
            self._read_from_misc = False
        try:
            # format-specific callback which will set headers and call
            #   add_synonym and read_taxon_line
            self.res.partition_parsing_fn(self)
            self._finish_partition_after_parse()
            self._partition_synonyms()
            self._populated = True
        except:
            self._read = False
            self._read_from_misc = None
            raise

    def flush(self):
        self.write_if_needed()

    def write_if_needed(self):
        if not self._populated
            return
        if self._subdirname_to_tp_roots:
            dh = self._misc_part
            dest = self.tax_fp_misc
            syndest = self.syn_fp_misc
            roots_file = os.path.join(self.tax_dir_misc, 'roots.txt')
        else:
            dh = self
            dest = self.tax_fp_unpartitioned
            syndest = self.syn_fp
            roots_file = os.path.join(self.tax_dir_unpartitioned, 'roots.txt')
        syn_id_order = _write_d_as_tsv(self.taxon_header, dh._id_to_line, dh._id_order, dest)
        if not dh._roots:
            _LOG.info('No root ids need to be written to "{}"'.format(roots_file))
        else:
            _LOG.info('Writing {} root_ids to "{}"'.format(len(dh._roots), roots_file))
            pd = os.path.split(roots_file)[0]
            assure_dir_exists(pd)
            with codecs.open(roots_file, 'w', encoding='utf-8') as outp:
                outp.write('\n'.join([str(i) for i in dh._roots]))
        if syndest is None:
            return
        _write_d_as_tsv(self.syn_header, dh._syn_by_id, syn_id_order, syndest)


def get_taxon_partition(res, fragment):
    ck = (TaxonPartition, res.id, fragment)
    c = TAX_SLICE_CACHE.get(ck)
    if c is not None:
        return c
    return TaxonPartition(res, fragment)

def _append_taxon(dict_to_write, id_order, dest_path):
    if not dict_to_write:
        _LOG.info('No records need to be appended to "{}"'.format(dest_path))
        return
    _LOG.info('Appending {} records to "{}"'.format(len(id_order), dest_path))
    ret = []
    with codecs.open(dest_path, 'a', encoding='utf-8') as outp:
        for i in id_order:
            el = dict_to_write.get(i)
            if el is not None:
                ret.append(i)
                outp.write(el)

        if len(id_order) != len(dict_to_write):
            oset = frozenset(ret)
            for key, line in dict_to_write.items():
                if key not in oset:
                    ret.append(key)
                    outp.write(line)
    return key

def _write_d_as_tsv(header, dict_to_write, id_order, dest_path):
    _LOG.info('Writing header to "{}"'.format(dest_path))
    pd = os.path.split(dest_path)[0]
    assure_dir_exists(pd)
    with codecs.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
    return _append_taxon(dict_to_write, id_order, dest_path)

'''
def _write_taxon_list(header, record_list, dest_path):
    _LOG.info('Writing header to "{}"'.format(dest_path))
    pd = os.path.split(dest_path)[0]
    assure_dir_exists(pd)
    with codecs.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
    _append_taxon_list(record_list, dest_path)


def _append_taxon_list(record_list, dest_path):
    if not record_list:
        _LOG.info('No records need to be appended to "{}"'.format(dest_path))
        return
    _LOG.info('Appending {} records to "{}"'.format(len(record_list), dest_path))
    with codecs.open(dest_path, 'a', encoding='utf-8') as outp:
        for line in record_list:
            outp.write(line)
'''