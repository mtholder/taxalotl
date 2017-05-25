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

    def try_del(self, ck):
        try:
            if ck in self._ck_to_obj:
                self.__delitem__(ck)
        except:
            _LOG.exception("caught and suppressed removal of key")
            if ck in self._ck_to_obj:
                del self._ck_to_obj[ck]

    def flush(self):
        kv = [(k, v) for k, v in self._ck_to_obj.items()]
        _ex = None
        for k, v in kv:
            if k in self._ck_to_obj:
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
        c = set(self._id_to_child_set.keys())
        c.update(self._id_to_line.keys())
        return c

    def _transfer_subtree(self, par_id, dest_part):
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
                dest_part.add_taxon(child_id, par_id, self._id_to_line[child_id])
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
        cids = set(self._id_to_child_set.keys())
        for dest_tax_part in other._root_to_lth.values():
            tpids = dest_tax_part.contained_ids()
            common = cids.intersection(tpids)
            if common:
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

    def _finish_partition_after_parse(self):
        """On entry _id_to_el will be set for the root elements (and some of their
            children), but taxa processed before their ancestors may have been missed.
        _id_to_child_list and _id_to_line are only filled for these
        """
        des_children_for_misc = []
        for uid, par_id in self._during_parse_root_to_par.items():
            line = self._id_to_line[uid]
            des_children_for_misc.append((uid, par_id, line))
        for uid, par_id in self._during_parse_root_to_par.items():
            match_el = self._root_to_lth[uid]
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
        #  _partition_synonyms that have now moved to the misc part
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
        if not self._populated:
            self._read_inputs()
        return self._misc_part.move_from_self_to_new_part(other)


class TaxonPartition(PartitionedTaxDirBase, PartitioningLightTaxHolder):
    def __init__(self, res, fragment):
        PartitioningLightTaxHolder.__init__(self, fragment)
        PartitionedTaxDirBase.__init__(self, res, fragment)
        self.treat_syn_as_taxa = self.syn_fp is None
        self._read = False
        self._read_from_partitioning_scratch = False
        self._read_from_misc = False
        self._fs_is_partitioned = None
        self._has_flushed = False

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
        '''
        Need to make this check for the input tax dir, not just scaffold...
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
        '''
        for subname, subroot in list_of_subdirname_and_roots:
            subfrag = os.path.join(self.fragment, subname)
            subtp = get_taxon_partition(self.res, subfrag)
            self._roots_for_sub.update(subroot)
            for r in subroot:
                self._root_to_lth[r] = subtp
            self._subdirname_to_tp_roots[subname] = (subtp, subroot)
        if self._populated:
            self._partition_from_in_mem()
        else:
            self._read_inputs()

    def _partition_from_in_mem(self):
        _LOG.info("_partition_from_in_mem for fragment \"{}\"".format(self.fragment))
        moved = set()
        assert not self._misc_part._populated
        for sub_tp, subroot in self._subdirname_to_tp_roots.values():
            _LOG.info("subroot {} for \"{}\"".format(subroot, sub_tp.fragment))
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

    def _move_data_to_empty_misc(self):
        assert not self._misc_part._populated
        m = self._misc_part
        for a in LightTaxonomyHolder._DATT:
            setattr(m, a, getattr(self, a))
            setattr(self, a, None)
        self._copy_shared_fields(m)
        m._populated = True

    def sub_tax_parts(self, include_misc=True):
        ret = [i for i in self._root_to_lth.values()]
        if include_misc:
            ret.append(self._misc_part)
        return ret

    def _copy_shared_fields(self, other):
        other.taxon_header = self.taxon_header
        other.syn_header = self.syn_header
        other.treat_syn_as_taxa = self.treat_syn_as_taxa


    def _read_inputs(self):
        self._read = True
        if os.path.exists(self.tax_fp_misc):
            self._read_from_partitioning_scratch = True
            self.tax_fp = self.tax_fp_misc
            self._read_from_misc = True
        else:
            self._read_from_partitioning_scratch = True
            self.tax_fp = self.tax_fp_unpartitioned
            self._read_from_misc = False
        try:
            # format-specific callback which will set headers and call
            #   add_synonym and read_taxon_line
            self.res.partition_parsing_fn(self)
            m = "prepart {} lines in misc. {} lines in self for {}"
            _LOG.info(m.format(len(self._misc_part._id_to_line),
                               len(self._id_to_line), self.fragment))
            self._finish_partition_after_parse()
            for el in self.sub_tax_parts():
                self._copy_shared_fields(el)
                el._populated = True
            self._populated = True
        except:
            self._read = False
            self._read_from_misc = None
            self._read_from_partitioning_scratch = False
            raise

    def flush(self):
        if self._has_flushed:
            return
        _LOG.info("flushing TaxonPartition for {}".format(self.fragment))
        wrote = self.write_if_needed()
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
            roots_file = os.path.join(self.tax_dir_misc, 'roots.txt')
        else:
            # _LOG.debug("write from self for {}".format(self.fragment))
            dh = self
            dest = self.tax_fp_unpartitioned
            syndest = self.syn_fp
            roots_file = os.path.join(self.tax_dir_unpartitioned, 'roots.txt')
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
