#!/usr/bin/env python
# from __future__ import print_function
import codecs
import os

from peyotl import get_logger, assure_dir_exists
from taxalotl.ott_schema import OTTTaxon

INP_TAXONOMY_DIRNAME = '__inputs__'
MISC_DIRNAME = '__misc__'
GEN_MAPPING_FILENAME = '__mapping__.json'
ROOTS_FILENAME = 'roots.txt'

_LOG = get_logger(__name__)

def get_root_ids_for_subset(tax_dir):
    rf = os.path.join(tax_dir, ROOTS_FILENAME)
    idset = set()
    if os.path.exists(rf):
        content = [int(i.strip()) for i in open(rf, 'r') if i.strip()]
        idset.update(content)
    return idset


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
        #_LOG.debug("TAX_SLICE_CACHE __delitem__ for {}".format(ck))
        obj = self._ck_to_obj.get(ck)
        if obj:
            del self._ck_to_obj[ck]
            #_LOG.debug("TAX_SLICE_CACHE call of flush for {}".format(ck))
            obj._flush()

    def try_del(self, ck):
        #_LOG.debug("TAX_SLICE_CACHE try_del for {}".format(ck))
        try:
            if ck in self._ck_to_obj:
                self.__delitem__(ck)
        except:
            _LOG.exception("caught and suppressed removal of key")
            if ck in self._ck_to_obj:
                del self._ck_to_obj[ck]

    def flush(self):
        #_LOG.debug("TAX_SLICE_CACHE flush")
        kv = [(k, v) for k, v in self._ck_to_obj.items()]
        self._ck_to_obj = {}
        _ex = None
        for k, v in kv:
            try:
                #_LOG.debug("TAX_SLICE_CACHE calling flush for {}".format(k))
                v._flush()
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

    def scaffold_tax_subdir_names(self):
        """Returns a list of subdirectory names for self.scaffold_dir with __misc__ suppressed"""
        if not os.path.isdir(self.scaffold_dir):
            return []
        n = []
        for x in os.listdir(self.scaffold_dir):
            if x == INP_TAXONOMY_DIRNAME or x == MISC_DIRNAME:
                continue
            if os.path.isdir(os.path.join(self.scaffold_dir, x)):
                n.append(x)
        return n


# noinspection PyProtectedMember
class LightTaxonomyHolder(object):
    _DATT = ['_id_order', '_id_to_line', '_id_to_child_set', '_syn_by_id', '_id_to_el', '_roots']

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
        self._populated = False
        self._has_unread_tax_inp = False

    def _del_data(self):
        for el in LightTaxonomyHolder._DATT:
            setattr(self, el, None)
        self._populated = False

    def add_taxon(self, uid, par_id, line):
        old = self._id_to_line.get(uid)
        assert old is None or old == line
        self._id_to_line[uid] = line
        self._id_order.append(uid)
        self._id_to_child_set.setdefault(par_id, set()).add(uid)

    add_moved_taxon = add_taxon
    add_taxon_from_higher_tax_part = add_taxon

    def contained_ids(self):
        c = set()
        if self._id_to_child_set:
            c.update(self._id_to_child_set.keys())
        if self._id_to_line:
            c.update(self._id_to_line.keys())
        return c

    def _transfer_subtree(self, par_id, dest_part):  # type (int, LightTaxonomyHolder) -> None
        child_set = self._id_to_child_set[par_id]
        # if self.fragment.startswith('Life/Archaea/Euryarchaeota'):
        #    _LOG.info("par_id {} -> child_list {}".format(par_id, child_set))
        self._id_to_el[par_id] = dest_part
        line = self._id_to_line.get(par_id)
        if line is not None:
            dest_part._id_to_line[par_id] = line
            del self._id_to_line[par_id]
        del self._id_to_child_set[par_id]
        dest_part._id_to_child_set.setdefault(par_id, set()).update(child_set)
        for child_id in child_set:
            self._id_to_el[child_id] = dest_part
            if child_id in self._id_to_child_set:
                self._transfer_subtree(child_id, dest_part)
            else:
                line = self._id_to_line.get(child_id)
                if line:
                    dest_part.add_taxon(child_id, par_id, line)
                    del self._id_to_line[child_id]

    def move_matched_synonyms(self, dest_tax_part):  # type: (PartitioningLightTaxHolder) -> None
        sk = set(self._syn_by_id.keys())
        sk.intersection_update(dest_tax_part.contained_ids())
        for s in sk:
            sd = self._syn_by_id[s]
            for pair in sd:
                dest_tax_part.add_synonym(s, pair[0], pair[1])
            del self._syn_by_id[s]

    def move_from_self_to_new_part(self, other):  # type: (PartitioningLightTaxHolder) -> None
        self._has_moved_taxa = True
        if other._has_unread_tax_inp:
            other._read_inputs(False)
        cids = set(self._id_to_child_set.keys())
        for dest_tax_part in other._root_to_lth.values():
            tpids = dest_tax_part.contained_ids()
            common = cids.intersection(tpids)
            if common:
                if dest_tax_part._has_unread_tax_inp:
                    dest_tax_part._read_inputs(False)
                for com_id in common:
                    if com_id in self._id_to_child_set:
                        m = "Transferring {} from {} to {}"
                        _LOG.info(m.format(com_id, self.fragment, dest_tax_part.fragment))
                        self._transfer_subtree(com_id, dest_tax_part)
                self.move_matched_synonyms(dest_tax_part)
                dest_tax_part._populated = True
                cids = set(self._id_to_child_set.keys())

    def add_synonym(self, accept_id, syn_id, line):
        if self.treat_syn_as_taxa:
            # CoL uses the taxonomy file for synonyms.
            assert syn_id is not None
            self.add_taxon(syn_id, None, line)
        else:
            self._syn_by_id.setdefault(accept_id, []).append((syn_id, line))


# noinspection PyProtectedMember
class PartitioningLightTaxHolder(LightTaxonomyHolder):
    def __init__(self, fragment):
        LightTaxonomyHolder.__init__(self, fragment)
        self._subdirname_to_tp_roots = {}
        self._misc_part = LightTaxonomyHolder(os.path.join(fragment, MISC_DIRNAME))
        self._roots_for_sub = set()
        self._root_to_lth = {}
        self._during_parse_root_to_par = {}
        self._has_moved_taxa = False # true when taxa have been moved to another partition

    def read_taxon_line(self, uid, par_id, line):
        if par_id:
            try:
                par_id = int(par_id)
            except:
                pass
        self._id_to_child_set.setdefault(par_id, set()).add(uid)
        self._id_to_line[uid] = line
        if uid in self._roots_for_sub:
            self._during_parse_root_to_par[uid] = par_id

    def sub_tax_parts(self, include_misc=True):
        ret = [i for i in self._root_to_lth.values()]
        if include_misc:
            # noinspection PyTypeChecker
            ret.append(self._misc_part)
        return ret

    def _finish_partition_after_parse(self):
        """On entry _id_to_el will be set for the root elements (and some of their
            children), but taxa processed before their ancestors may have been missed.
        _id_to_child_list and _id_to_line are only filled for these
        """
        _LOG.debug("_finish_partition_after_parse for {}".format(self.fragment))
        for tp in self.sub_tax_parts(include_misc=False):
            if tp._has_unread_tax_inp:
                tp._read_inputs(False)
        des_children_for_misc = []
        for uid, par_id in self._during_parse_root_to_par.items():
            line = self._id_to_line[uid]
            des_children_for_misc.append((uid, par_id, line))
        for uid, par_id in self._during_parse_root_to_par.items():
            match_el = self._root_to_lth[uid]
            m = 'transferring taxon {} from "{}" to "{}"'
            _LOG.debug(m.format(uid, self.fragment, match_el.fragment))
            match_el._roots.add(uid)
            if uid in self._id_to_child_set:
                self._transfer_subtree(uid, match_el)
            elif uid in self._id_to_line:
                match_el._id_to_line[uid] = self._id_to_line[uid]
                del self._id_to_line[uid]
            pc = self._id_to_child_set.get(par_id)
            if pc:
                pc.remove(uid)
            self._id_to_el[uid] = match_el
        assert not self._misc_part._id_to_child_set
        assert not self._misc_part._id_to_line
        assert not self._misc_part._id_to_el
        # Move all data to misc, but make a copy of th id to el that we'l
        self._move_data_to_empty_misc()
        for el in des_children_for_misc:
            self._misc_part.add_moved_taxon(el[0], el[1], el[2])
        # _partition_synonyms that have now moved to the misc part
        to_del = set()
        for accept_id, i_l_list in self._misc_part._syn_by_id.items():
            match_el = self._misc_part._id_to_el.get(accept_id)
            if match_el is not None:
                for syn_id, line in i_l_list:
                    match_el.add_synonym(accept_id, syn_id, line)
                to_del.add(accept_id)
        for i in to_del:
            del self._misc_part._syn_by_id[i]

    def _read_inputs(self):
        raise NotImplementedError("_read_input pure virtual in PartitioningLightTaxHolder")

    def move_from_misc_to_new_part(self, other):  # type: (PartitioningLightTaxHolder) -> None
        self._has_moved_taxa = True
        if not self._populated:
            self._read_inputs()
        return self._misc_part.move_from_self_to_new_part(other)

    def _move_data_to_empty_misc(self):
        assert not self._misc_part._populated
        m = self._misc_part
        for a in LightTaxonomyHolder._DATT:
            setattr(m, a, getattr(self, a))
            setattr(self, a, None)
        self._copy_shared_fields(m)
        m._populated = True

    def _copy_shared_fields(self, other):
        other.taxon_header = self.taxon_header
        other.syn_header = self.syn_header
        other.treat_syn_as_taxa = self.treat_syn_as_taxa


# noinspection PyProtectedMember
class TaxonPartition(PartitionedTaxDirBase, PartitioningLightTaxHolder):
    def __init__(self, res, fragment):
        PartitioningLightTaxHolder.__init__(self, fragment)
        PartitionedTaxDirBase.__init__(self, res, fragment)
        self.treat_syn_as_taxa = self.syn_fp is None
        self._read_from_fs = False
        self._read_from_partitioning_scratch = False
        self._read_from_misc = False
        self._fs_is_partitioned = None
        self._has_flushed = False
        self._external_inp_fp = None

    @property
    def external_input_fp(self):
        return self._external_inp_fp

    @external_input_fp.setter
    def external_input_fp(self, value):
        self._external_inp_fp = value

    def _diagnose_state_of_fs(self):
        if os.path.exists(self.tax_fp_misc):
            self._fs_is_partitioned = True
        elif os.path.exists(self.tax_fp_unpartitioned):
            self._fs_is_partitioned = False
        else:
            self._fs_is_unpartitioned = None

    def taxa_files_exist_for_a_frag(self, frag):
        if os.path.exists(self.res.get_taxon_filepath_for_part(frag)):
            return True
        return os.path.exists(self.res.get_misc_taxon_dir_for_part(frag))

    def do_partition(self, list_of_subdirname_and_roots):
        if self._subdirname_to_tp_roots:
            raise ValueError("do_partition called twice for {}".format(self.fragment))
        if not self._populated:
            self._diagnose_state_of_fs()
            if self._fs_is_partitioned is None:
                m = "Taxa files not found for {} and TaxonPartition is empty"
                _LOG.info(m.format(self.fragment))
        cur_sub_names = self.scaffold_tax_subdir_names()
        do_part_if_reading = True
        having_inp_to_read = set()
        if cur_sub_names:
            req_fulfilled = True
            some_part_found = False
            for x in list_of_subdirname_and_roots:
                subname = x[0]
                if subname in cur_sub_names:
                    subfrag = os.path.join(self.fragment, subname)
                    if self.taxa_files_exist_for_a_frag(subfrag):
                        _LOG.info("previous content for {}".format(subfrag))
                        some_part_found = True
                        having_inp_to_read.add(subname)
                    else:
                        _LOG.warn("no previous taxonomic content for {}".format(subfrag))
                        req_fulfilled = False
                else:
                    _LOG.warn("no previous subdir for {}".format(subname))
                    req_fulfilled = False
            if some_part_found:
                do_part_if_reading = False
                quant = 'All' if req_fulfilled else 'Some'
                m = "{} subdir partitions found for {}. No more partitioning of parent will be done!"
                _LOG.warn(m.format(quant, self.fragment))
        for subname, subroot in list_of_subdirname_and_roots:
            subfrag = os.path.join(self.fragment, subname)
            subtp = get_taxon_partition(self.res, subfrag)
            if subname in having_inp_to_read:
                subtp._has_unread_tax_inp = True
            self._roots_for_sub.update(subroot)
            for r in subroot:
                self._root_to_lth[r] = subtp
            self._subdirname_to_tp_roots[subname] = (subtp, subroot)
        if self._populated:
            self._partition_from_in_mem()
        else:
            self._read_inputs(do_part_if_reading)

    def _partition_from_in_mem(self):
        _LOG.info("_partition_from_in_mem for fragment \"{}\"".format(self.fragment))
        if self._misc_part._populated:
            m = "_partition_from_in_mem called for {}, but misc already has {}"
            raise ValueError(m.format(self.fragment, self.contained_ids()))
        self._has_moved_taxa = True
        for sub_tp, subroot in self._subdirname_to_tp_roots.values():
            _LOG.info("subroot {} for \"{}\"".format(subroot, sub_tp.fragment))
            if sub_tp._has_unread_tax_inp:
                sub_tp._read_inputs(False)
            x = self._id_to_child_set
            # if self.fragment.startswith('Life/Archaea/Euryarchaeota'):
            #    _LOG.info(" self._id_to_child_set = {}".format(repr(x)))
            #    #_LOG.info(" self._id_to_line = {}".format(self._id_to_line))
            for r in subroot:
                # if self.fragment.startswith('Life/Archaea/Euryarchaeota'):
                #     _LOG.info(" checking subroot {}".format(repr(r)))
                if r in x:
                    # if self.fragment.startswith('Life/Archaea/Euryarchaeota'):
                    #     _LOG.info(" {} in _id_to_child_set".format(repr(r)))
                    self._transfer_subtree(r, sub_tp)
                    sub_tp._roots.add(r)
                elif r in self._id_to_line:
                    # if self.fragment.startswith('Life/Archaea/Euryarchaeota'):
                    #     _LOG.info(" {} in _id_to_line".format(repr(r)))
                    sub_tp._id_to_line[r] = self._id_to_line[r]
                    sub_tp._roots.add(r)
                    # else:
                    #     if self.fragment.startswith('Life/Archaea/Euryarchaeota'):
                    #         _LOG.info(" not stored".format(r))

            self.move_matched_synonyms(sub_tp)
            self._copy_shared_fields(sub_tp)
            sub_tp._populated = True
        self._move_data_to_empty_misc()

    def read_inputs_for_read_only(self):
        # Only to be used for accessors
        assert not self._read_from_fs
        assert not self._populated
        self._read_inputs(do_part_if_reading=False)

    def get_root_ids(self):
        return set(self._roots)

    def get_id_to_ott_taxon(self):
        id_to_obj = {}
        for line in self._id_to_line.values():
            obj = OTTTaxon(line)
            oid = obj.id
            assert oid not in id_to_obj
            id_to_obj[oid] = obj
        return id_to_obj

    def _read_inputs(self, do_part_if_reading=True):
        self._read_from_fs = True
        self._has_unread_tax_inp = False

        if self._external_inp_fp:
            self._read_from_partitioning_scratch = True
            self.tax_fp = self._external_inp_fp
            tax_dir = os.path.split(self.tax_fp)[0]
        else:
            if os.path.exists(self.tax_fp_misc):
                self._read_from_partitioning_scratch = True
                self.tax_fp = self.tax_fp_misc
                tax_dir = self.tax_dir_misc
                self._read_from_misc = True
            else:
                self._read_from_partitioning_scratch = True
                self.tax_fp = self.tax_fp_unpartitioned
                self._read_from_misc = False
                tax_dir = self.tax_dir_unpartitioned
        try:
            # format-specific callback which will set headers and call
            #   add_synonym and read_taxon_line
            _LOG.debug("About to parse taxa. Looking for root ids: {}".format(self._roots_for_sub))
            self.res.partition_parsing_fn(self)
            read_roots = get_root_ids_for_subset(os.path.join(tax_dir, ROOTS_FILENAME))
            self._roots.update(read_roots)
            m = "prepart {} taxa in {}"
            _LOG.info(m.format(len(self._misc_part._id_to_line) + len(self._id_to_line), self.fragment))
            self._read_from_fs
            if do_part_if_reading:
                self._has_moved_taxa = True
                self._finish_partition_after_parse()
                for el in self.sub_tax_parts():
                    self._copy_shared_fields(el)
                    el._populated = True
            self._populated = True
        except:
            self._read_from_fs = False
            self._read_from_misc = None
            self._read_from_partitioning_scratch = False
            raise

    def _flush(self):
        if self._has_flushed:
            _LOG.info("duplicate flush of TaxonPartition for {} ignored.".format(self.fragment))
            return
        if not self._has_moved_taxa:
            if self._read_from_fs:
                _LOG.info("Flush of unaltered TaxonPartition for {} ignored".format(self.fragment))
                return
        _LOG.info("flushing TaxonPartition for {}".format(self.fragment))
        self.write_if_needed()
        if self._read_from_misc is False and self._read_from_partitioning_scratch:
            tr = [self.tax_fp_unpartitioned]
            if self.syn_fp:
                tr.append(self.syn_fp)
            for f in tr:
                if os.path.exists(f):
                    _LOG.info("removing pre-partitioned file at {}".format(f))
                    try:
                        os.unlink(f)
                    except:
                        _LOG.exception("could not remove {}".format(f))
        self._has_flushed = True
        TAX_SLICE_CACHE.try_del(self.cache_key)
        self._del_data()

    def write_if_needed(self):
        if not self._populated:
            _LOG.info("write not needed for {} not populated".format(self.fragment))
            return False
        if self._subdirname_to_tp_roots:
            # _LOG.debug("write from misc for {}".format(self.fragment))
            dh = self._misc_part
            dest = self.tax_fp_misc
            syndest = self.syn_fp_misc
            roots_file = os.path.join(self.tax_dir_misc, ROOTS_FILENAME)
        else:
            # _LOG.debug("write from self for {}".format(self.fragment))
            dh = self
            dest = self.tax_fp_unpartitioned
            syndest = self.syn_fp
            roots_file = os.path.join(self.tax_dir_unpartitioned, ROOTS_FILENAME)
        if not dh._id_to_line:
            _LOG.debug("write not needed for {} no records".format(self.fragment))
            return False
        syn_id_order = _write_d_as_tsv(self.taxon_header, dh._id_to_line, dh._id_order, dest)
        if not dh._roots:
            _LOG.debug('No root ids need to be written to "{}"'.format(roots_file))
        else:
            _LOG.debug('Writing {} root_ids to "{}"'.format(len(dh._roots), roots_file))
            pd = os.path.split(roots_file)[0]
            assure_dir_exists(pd)
            with codecs.open(roots_file, 'w', encoding='utf-8') as outp:
                outp.write('\n'.join([str(i) for i in dh._roots]))
        if syndest is None:
            return
        _write_syn_d_as_tsv(self.syn_header, dh._syn_by_id, syn_id_order, syndest)
        return True


def get_taxon_partition(res, fragment):
    ck = (TaxonPartition, res.id, fragment)
    c = TAX_SLICE_CACHE.get(ck)
    if c is not None:
        return c
    return TaxonPartition(res, fragment)


def _write_d_as_tsv(header, dict_to_write, id_order, dest_path):
    if not dict_to_write:
        return
    ret = []
    pd = os.path.split(dest_path)[0]
    assure_dir_exists(pd)
    _LOG.info('Writing {} records to "{}"'.format(len(dict_to_write), dest_path))
    with codecs.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
        for i in id_order:
            el = dict_to_write.get(i)
            if el is not None:
                ret.append(i)
                outp.write(el)
        oset = frozenset(ret)
        for key, line in dict_to_write.items():
            if key not in oset:
                ret.append(key)
                outp.write(line)
    return ret


def _write_syn_d_as_tsv(header, dict_to_write, id_order, dest_path):
    ltw = []
    for i in id_order:
        synlist = dict_to_write.get(i)
        if synlist is not None:
            for p in synlist:
                ltw.append(p[1])
    oset = frozenset(id_order)
    for key, synlist in dict_to_write.items():
        if key not in oset:
            for syn_pair in synlist:
                ltw.append(syn_pair[1])
    if not ltw:
        return
    x = len(ltw)
    pd = os.path.split(dest_path)[0]
    assure_dir_exists(pd)
    _LOG.info('Writing {} records to "{}"'.format(x, dest_path))
    with codecs.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
        for l in ltw:
            outp.write(l)


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
