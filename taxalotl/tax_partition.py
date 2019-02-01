#!/usr/bin/env python
import io
import os
from copy import copy
from contextlib import contextmanager

from peyotl import get_logger, assure_dir_exists, read_as_json, write_as_json

from .taxon import Taxon
from .tree import TaxonForest
from .ott_schema import HEADER_TO_LINE_PARSER
from .util import unlink

INP_TAXONOMY_DIRNAME = '__inputs__'
MISC_DIRNAME = '__misc__'
GEN_MAPPING_FILENAME = '__mapping__.json'
ROOTS_FILENAME = '__roots__.json'
ACCUM_DES_FILENAME = '__accum_des__.json'

_LOG = get_logger(__name__)


def get_taxonomies_for_dir(taxon_dir):
    inps_dir = os.path.join(taxon_dir, INP_TAXONOMY_DIRNAME)
    resource_ids = []
    for name in os.listdir(inps_dir):
        fp = os.path.join(inps_dir, name)
        if os.path.isdir(fp) and os.path.isfile(os.path.join(fp, ROOTS_FILENAME)):
            resource_ids.append(name)
    return resource_ids


def get_roots_for_subset(tax_dir, misc_tax_dir):
    return _read_json_and_coerce_to_otttaxon(tax_dir, misc_tax_dir, ROOTS_FILENAME)


def get_accum_des_for_subset(tax_dir, misc_tax_dir):
    return _read_json_and_coerce_to_otttaxon(tax_dir, misc_tax_dir, ACCUM_DES_FILENAME)


def _read_json_and_coerce_to_otttaxon(tax_dir, misc_tax_dir, fn):
    r = {}
    for td in [tax_dir, misc_tax_dir]:
        rf = os.path.join(td, fn)
        if os.path.exists(rf):
            rd = read_as_json(rf)
            for k, v in rd.items():
                try:
                    k = int(k)
                except:
                    pass
                r[k] = Taxon(d=v)
    return r


# noinspection PyProtectedMember
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
            obj._flush()

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
        self._ck_to_obj = {}
        _ex = None
        for k, v in kv:
            try:
                v._flush()
            except Exception as x:
                _LOG.exception('exception in flushing')
                _ex = x
        if _ex is not None:
            raise _ex

    # noinspection PyMethodMayBeStatic
    def get_taxon_partition(self, res, fragment):
        return get_taxon_partition(res, fragment)

    def clear_without_flush(self, ck):
        if ck in self._ck_to_obj:
            del self._ck_to_obj[ck]


TAX_SLICE_CACHE = TaxonomySliceCache()


@contextmanager
def use_tax_partitions():
    yield TAX_SLICE_CACHE
    TAX_SLICE_CACHE.flush()


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
        self.synonyms_filename = self.res.synonyms_filename
        self.cache_key = (self.__class__, self.src_id, self.fragment)
        assert TAX_SLICE_CACHE.get(self.cache_key) is None
        TAX_SLICE_CACHE[self.cache_key] = self

    @property
    def input_taxdir(self):
        return os.path.split(self.tax_fp)[0]

    @property
    def input_synonyms_filepath(self):
        if self.synonyms_filename:
            return os.path.join(self.input_taxdir, self.synonyms_filename)
        return None

    @property
    def output_synonyms_filepath(self):
        if not self.synonyms_filename:
            return None
        pd = self.tax_dir_misc if self._subdirname_to_tp_roots else self.tax_dir_unpartitioned
        return os.path.join(pd, self.synonyms_filename)

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


class Synonym(object):
    def __init__(self, valid_tax_id, name, syn_type=None, syn_id=None):
        if syn_type is None:
            syn_type = 'synonym'
        self.valid_tax_id = valid_tax_id
        self.name = name
        self.syn_type = syn_type
        self.syn_id = syn_id

    def to_serializable_dict(self):
        d = {'valid_tax_id': self.valid_tax_id,
             'name': self.name
             }
        if self.syn_id:
            d['synonym_id'] = self.syn_id
        if self.syn_type != 'synonym':
            d['synonym'] = self.syn_type
        return d

    def __repr__(self):
        sis = ', syn_id={}'.format(self.syn_id) if self.syn_id else ''
        tis = ', syn_type={}'.format(repr(self.syn_type)) if self.syn_type != 'synonym' else ''
        m = 'Synonym({}, {}{}{})'
        return m.format(self.valid_tax_id, repr(self.name), tis, sis)

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        return hash((self.syn_id, self.valid_tax_id, self.name, self.syn_type))

    def __eq__(self, other):
        if isinstance(other, Synonym):
            return self.valid_tax_id == other.valid_tax_id \
                   and self.name == other.name \
                   and self.syn_type == other.syn_type \
                   and self.syn_id == self.syn_id


_VALID_SYN_TYPES = {'acronym', 'ambiguous synonym',
                    'authority',
                    'blast name',
                    'common name',
                    'equivalent name',
                    'genbank common name',
                    'genbank synonym',
                    'homotypic synonym',
                    'includes', 'misapplied name',
                    'misnomer',
                    'misspelling',
                    'orthographia',
                    'proparte synonym',
                    'synonym',
                    'type material',
                    }
IGNORE_SYN_TYPES = {'acronym',
                    'authority',
                    'blast name',
                    'common name',
                    'genbank common name',
                    'genbank synonym',
                    'type material',
                    }


class SynonymInterpreter(object):
    def __init__(self, header):
        if header.endswith('\n'):
            header = header[:-1]
        self.fields = [i.strip() for i in header.split('\t|\t') if i.strip()]
        assert 'uid' in self.fields
        self._uid_ind = self.fields.index('uid')
        self._name_ind = self.fields.index('name')
        self._type_ind = self.fields.index('type')

    def interpret(self, uid, syn_id_line_tuple):
        syn_id, line = syn_id_line_tuple
        sl = line.split('\t|\t')
        suid = sl[self._uid_ind]
        name = sl[self._name_ind].strip()
        syn_type = sl[self._type_ind].strip()
        if syn_type not in _VALID_SYN_TYPES:
            m = 'synonym_type "{}" not recognized in for ({}, "{}")'
            raise ValueError(m.format(syn_type, uid, name))
        assert uid == int(suid)
        return Synonym(valid_tax_id=uid, name=name, syn_type=syn_type, syn_id=syn_id)


# noinspection PyProtectedMember
class LightTaxonomyHolder(object):
    _DATT = ['_des_in_other_slices',
             '_id_order',
             '_id_to_child_set',
             '_id_to_el',
             '_id_to_line',
             '_syn_by_id',
             '_roots',
             ]

    def __init__(self, fragment):
        self.fragment = fragment
        self._id_order = []
        self._id_to_line = {}  # id -> line
        self._id_to_child_set = {}  # id -> set of child IDs
        self._id_to_el = {}
        self._roots = {}
        self._des_in_other_slices = {}
        self._syn_by_id = {}  # accepted_id -> list of synonym lines
        self._parsed_syn_by_id = None
        self.taxon_header = None
        self.syn_header = None
        self.treat_syn_as_taxa = False
        self._populated = False
        self._has_unread_tax_inp = False
        self._has_moved_taxa = False  # true when taxa have been moved to another partition

    @property
    def write_taxon_header(self):
        from taxalotl.ott_schema import INP_FLAGGED_OTT_TAXONOMY_HEADER
        return INP_FLAGGED_OTT_TAXONOMY_HEADER

    @property
    def synonyms_by_id(self):
        return copy(self._syn_by_id)

    def parsed_synonyms_by_id(self, ignored_syn_types=None):
        if self._parsed_syn_by_id is None:
            p = {}
            if self._syn_by_id:
                si = SynonymInterpreter(self.syn_header)
                for uid, line_stub_list in self._syn_by_id.items():
                    ps = set()
                    for i in line_stub_list:
                        syn = si.interpret(uid, i)
                        st = syn.syn_type
                        if ignored_syn_types is None or st not in ignored_syn_types:
                            ps.add(syn)
                    if ps:
                        p[uid] = ps
            self._parsed_syn_by_id = p
        return copy(self._parsed_syn_by_id)

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

    add_taxon_from_higher_tax_part = add_taxon

    def contained_ids(self):
        c = set()
        if self._id_to_child_set:
            c.update(self._id_to_child_set.keys())
        if self._id_to_line:
            c.update(self._id_to_line.keys())
        return c

    def _add_root(self, uid, taxon):
        self._roots[uid] = taxon.to_serializable_dict()

    def line_to_taxon(self, line=None, uid=None):
        if line is None:
            line = self._id_to_line[uid]
        return Taxon(line, line_parser=HEADER_TO_LINE_PARSER[self.taxon_header])

    def _transfer_line(self, uid, dest_part,
                       as_root=False):  # type (int, LightTaxonomyHolder, bool) -> None
        line = self._id_to_line[uid]
        taxon = self.line_to_taxon(line)
        if as_root:
            dest_part._add_root(uid, taxon)
        d = taxon.to_serializable_dict()
        d['fragment'] = dest_part.fragment
        self._des_in_other_slices[uid] = d
        dest_part._id_to_line[uid] = line
        del self._id_to_line[uid]

    def _transfer_subtree(self, par_id, dest_part,
                          as_root=False):  # type (int, LightTaxonomyHolder) -> None
        self._has_moved_taxa = True
        taxon = self.line_to_taxon(uid=par_id)
        if as_root:
            dest_part._add_root(par_id, taxon)
        d = taxon.to_serializable_dict()
        d['fragment'] = dest_part.fragment
        self._des_in_other_slices[par_id] = d
        self._transfer_subtree_rec(par_id, dest_part)

    def _transfer_subtree_rec(self, par_id, dest_part):  # type (int, LightTaxonomyHolder) -> None
        assert self is not dest_part
        assert self.fragment != dest_part.fragment
        child_set = self._id_to_child_set[par_id]
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
                self._transfer_subtree_rec(child_id, dest_part)
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
        lth_frag_to_root_id_set = {}
        for root_id, dest_tax_part in other._root_to_lth.items():
            lth_frag_to_root_id_set.setdefault(dest_tax_part.fragment, set()).add(root_id)

        for dest_tax_part in other._root_to_lth.values():
            tpids = set(dest_tax_part.contained_ids())
            tpids.update(lth_frag_to_root_id_set[dest_tax_part.fragment])
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
        ls = fragment.split('/')
        if len(ls) > 1:
            assert ls[-2] != ls[-1]
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
        if uid in self._id_to_line:
            raise ValueError("Repeated uid {} in line {}".format(uid, line))
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
        for tp in self.sub_tax_parts(include_misc=False):
            if tp._has_unread_tax_inp:
                tp._read_inputs(False)
        for uid, par_id in self._during_parse_root_to_par.items():
            match_el = self._root_to_lth[uid]
            if uid in self._id_to_child_set:
                self._transfer_subtree(uid, match_el, as_root=True)
            elif uid in self._id_to_line:
                self._transfer_line(uid, match_el, as_root=True)
            pc = self._id_to_child_set.get(par_id)
            if pc:
                pc.remove(uid)
            self._id_to_el[uid] = match_el
        assert not self._misc_part._id_to_child_set
        assert not self._misc_part._id_to_line
        assert not self._misc_part._id_to_el
        # Move all data to misc, but make a copy of th id to el that we'l
        self._move_data_to_empty_misc()
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

    def _read_inputs(self, do_part_if_reading=True):
        raise NotImplementedError("_read_input pure virtual in PartitioningLightTaxHolder")

    def move_from_misc_to_new_part(self, other):
        self._has_moved_taxa = True
        if not self._populated:
            self._read_inputs()
        if not self._misc_part._populated:
            self._move_data_to_empty_misc()
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
        self.treat_syn_as_taxa = self.synonyms_filename is None
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
            if (self._fs_is_partitioned is None) and (not self._external_inp_fp):
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
                m = "{} subdir partitions found for {}. No more partitioning will be done!"
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
            if sub_tp._has_unread_tax_inp:
                sub_tp._read_inputs(False)
            x = self._id_to_child_set
            for r in subroot:
                if r in x:
                    self._transfer_subtree(r, sub_tp, as_root=True)
                elif r in self._id_to_line:
                    self._transfer_line(r, sub_tp, as_root=True)
                else:
                    pass
            self.move_matched_synonyms(sub_tp)
            self._copy_shared_fields(sub_tp)
            sub_tp._populated = True
        self._move_data_to_empty_misc()

    def _debug_validity_check(self):
        self.read_inputs_for_read_only()
        _LOG.debug('{} roots = {}'.format(self.fragment, self._roots))
        id_to_par = {}
        errs = []
        warnings = []
        for par_id, cs in self._id_to_child_set.items():
            for c in cs:
                id_to_par[c] = par_id
        known_id_set = set(self._id_to_line.keys())
        roots_set = set(self._roots.keys())
        for uid in self._id_to_line.keys():
            par_id = id_to_par.get(uid)
            if not par_id:
                if uid not in self._roots:
                    m = 'ID {} does not have a parent in this slice, but is not listed in the roots'
                    if self.fragment == 'Life':
                        warnings.append(m.format(uid))
                    else:
                        errs.append(m.format(uid))
                    roots_set.add(uid)
        _LOG.debug('{} elements in self._id_to_line'.format(len(self._id_to_line)))
        for uid in self._syn_by_id.keys():
            if uid not in self._id_to_line:
                m = 'synonyms ID for {}, in syn_by_id but not in id_to_line'.format(uid)
                errs.append(m)
        known_id_set.update(self._roots.keys())
        for k in self._roots.keys():
            if k not in self._id_to_line:
                m = 'root ID {}, but roots not in id_to_line'.format(k)
                errs.append(m)
        for k, v in self._des_in_other_slices.items():
            if k in self._id_to_line:
                m = 'ID {} flagged as being in another slice, but it is still in id_to_line'.format(
                    k)
                errs.append(m)
            if v.get('par_id'):
                if v.get('par_id') not in self._id_to_line:
                    m = 'slice does not hold parent {} of id {} which is flagged as being in another slice'
                    m = m.format(v.get('par_id'), k)
                    warnings.append(m)
        if warnings:
            m = '{} warning(s): {}'.format(len(warnings), '\n'.join(warnings))
            _LOG.warn(m)
        if errs:
            m = '{} error(s): {}'.format(len(errs), '\n'.join(errs))
            raise ValueError(m)
        return roots_set, known_id_set

    def _debug_check_subtree_ids(self, root_id_set, all_id_set):
        self.read_inputs_for_read_only()
        self_ids_set = set()
        to_deal_with = set(root_id_set)
        unrecognized_set = set()
        while to_deal_with:
            ntdw = set()
            for r in to_deal_with:
                cs = self._id_to_child_set.get(r)
                if cs is None:
                    if r not in self._id_to_line:
                        unrecognized_set.add(r)
                elif cs:
                    self_ids_set.add(r)
                    self_ids_set.update(cs)
                    ntdw.update(cs)
            to_deal_with = ntdw
        missed_ids = self_ids_set - all_id_set
        extra_ids = all_id_set - self_ids_set
        errs = []
        # _LOG.debug('root_id_set = {}'.format(root_id_set))
        # _LOG.debug('len(self._id_to_child_set) = {} len(self_ids_set) = {} '.format(
        #     len(self._id_to_child_set), len(self_ids_set)))
        if unrecognized_set:
            x = list(unrecognized_set)
            x.sort()
            mind = len(x) if len(x) < 100 else 100
            m = 'IDs not known to {} read from {}: {}'
            errs.append(m.format(self.fragment, self.tax_fp, x[:mind]))
        if missed_ids:
            x = list(missed_ids)
            x.sort()
            mind = len(x) if len(x) < 100 else 100
            m = 'IDs expected in subtree according to {} read from {}, but not found: {}'
            errs.append(m.format(self.fragment, self.tax_fp, x[:mind]))
        if extra_ids:
            x = list(extra_ids)
            x.sort()
            mind = len(x) if len(x) < 100 else 100
            m = 'IDs included in subtree, but not expected by {} read from {}: {}'
            errs.append(m.format(self.fragment, self.tax_fp, x[:mind]))
        if errs:
            m = '{} error(s): {}'.format(len(errs), '\n'.join(errs))
            raise ValueError(m)

    def read_inputs_for_read_only(self):
        # Only to be used for accessors
        if not self._read_from_fs:
            assert not self._populated
            self._read_inputs(do_part_if_reading=False)

    def get_root_ids(self):
        return set(self._roots.keys())

    def get_id_to_ott_taxon(self):
        id_to_obj = {}
        lp = HEADER_TO_LINE_PARSER[self.taxon_header]
        for line in self._id_to_line.values():
            obj = Taxon(line, line_parser=lp)
            oid = obj.id
            assert oid not in id_to_obj
            id_to_obj[oid] = obj
        return id_to_obj

    def get_taxa_as_forest(self):
        return TaxonForest(id_to_taxon=self.get_id_to_ott_taxon(), taxon_partition=self)

    def active_tax_dir(self):
        if self._populated:
            return os.path.split(self.tax_fp)[0]
        raise NotImplementedError('active_tax_dir on unpopulated')

    def _read_inputs(self, do_part_if_reading=True):
        self._has_unread_tax_inp = False

        if self._external_inp_fp:
            self._read_from_partitioning_scratch = True
            self.tax_fp = self._external_inp_fp
        else:
            if os.path.exists(self.tax_fp_misc):
                self._read_from_partitioning_scratch = True
                self.tax_fp = self.tax_fp_misc
                self._read_from_misc = True
            else:
                self._read_from_partitioning_scratch = True
                self.tax_fp = self.tax_fp_unpartitioned
                self._read_from_misc = False
        try:
            self.res.partition_parsing_fn(self)
            read_roots = self._read_roots()
            self._roots.update(read_roots)
            self._des_in_other_slices.update(self.read_acccumulated_des())
            self._read_from_fs = True
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

    def _read_roots(self):
        return get_roots_for_subset(self.tax_dir_unpartitioned, self.tax_dir_misc)

    def read_acccumulated_des(self):
        return get_accum_des_for_subset(self.tax_dir_unpartitioned, self.tax_dir_misc)

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
            if self.output_synonyms_filepath:
                tr.append(self.output_synonyms_filepath)
            tr.append(os.path.join(self.tax_dir_unpartitioned, ACCUM_DES_FILENAME))
            for f in tr:
                if os.path.exists(f):
                    try:
                        unlink(f)
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
            out_dir = self.tax_dir_misc
        else:
            # _LOG.debug("write from self for {}".format(self.fragment))
            dh = self
            dest = self.tax_fp_unpartitioned
            out_dir = self.tax_dir_unpartitioned
        roots_file = os.path.join(out_dir, ROOTS_FILENAME)
        if not dh._id_to_line:
            _LOG.debug("write not needed for {} no records".format(self.fragment))
            syn_id_order = []
        else:
            syn_id_order = _write_d_as_tsv(self.write_taxon_header, dh._id_to_line, dh._id_order,
                                           dest)
        if not dh._roots:
            _LOG.debug('No root ids need to be written to "{}"'.format(roots_file))
        else:
            _LOG.debug('Writing {} root_ids to "{}"'.format(len(dh._roots), roots_file))
            write_taxon_json(dh._roots, roots_file)
        syndest = self.output_synonyms_filepath
        if syndest is not None:
            _write_syn_d_as_tsv(self.syn_header, dh._syn_by_id, syn_id_order, syndest)
        if dh._des_in_other_slices:
            write_taxon_json(dh._des_in_other_slices, os.path.join(out_dir, ACCUM_DES_FILENAME))
        return True


def write_taxon_json(obj, filepath):
    out_dir = os.path.split(filepath)[0]
    if out_dir:
        assure_dir_exists(out_dir)
    dtw = {}
    for k, v in obj.items():
        if isinstance(v, Taxon):
            dtw[k] = v.to_serializable_dict()
        else:
            dtw[k] = v
    write_as_json(dtw, filepath, separators=(',', ": "), indent=1)


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
    with io.open(dest_path, 'w', encoding='utf-8') as outp:
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
    if dict_to_write:
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
    with io.open(dest_path, 'w', encoding='utf-8') as outp:
        outp.write(header)
        for l in ltw:
            outp.write(l)
