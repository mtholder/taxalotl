#!/usr/bin/env python3

from sqlalchemy.orm import sessionmaker
import re
from peyotl import get_logger
from datetime import datetime

from .col_db_model import (AccordingTo, Authority,
                           CoLLink,
                           Dataset, Description,
                           Edge,
                           Modified,
                           NameFragment,
                           Rank,
                           SciNameUID, Synonym, SynonymStatus,
                           Taxon, TaxonName, TaxonConcept, TaxonSuppInfo, TaxUUID,
                           ValidNameStatus,
                           )

from ..taxonomic_ranks import (is_higher_taxon, is_infraspecific, RANK_TO_SORTING_NUMBER)

_LOG = get_logger(__name__)


def replace_values_with_keys(engine, parser, synonyms, by_par, par_less):
    def gen_iter():
        return gen_anc_first_then_syn(parser, synonyms, by_par, par_less)

    _status_update('assigning traversal indices')
    assign_traversals(parser, by_par, par_less)
    _status_update('dataset')
    gen_dataset_fk(engine, parser, gen_iter())
    _status_update('modtime')
    gen_modified_fk(engine, parser, gen_iter())
    _status_update('according_to')
    gen_according_to_fk(engine, parser, gen_iter())
    _status_update('description')
    gen_description_fk(engine, parser, gen_iter())
    _status_update('taxon_concept')
    gen_taxon_concept_fk(engine, parser, gen_iter())
    _status_update('sci_name_id')
    gen_sci_name_fk(engine, parser, gen_iter())
    _status_update('uuid and ColLink')
    gen_uuid_and_col_link_fk(engine, parser, gen_iter())
    _status_update('authority')
    gen_authority_fk(engine, parser, gen_iter())
    _status_update('name_status')
    gen_name_status_fk(engine, parser, gen_iter())
    _status_update('taxon_supplemental_info')
    append_tax_sup_inf_fk(engine, parser, gen_iter())
    _status_update('taxon_name')
    gen_name_fk(engine, parser, gen_iter())
    _status_update('taxon_rank')
    gen_rank_fk(engine, parser, gen_anc_first_for_each_subtree(parser, by_par, par_less))
    _status_update('taxon')
    tax_id2_tax_pk = gen_taxon_fk(engine, parser,
                                  gen_anc_first_then_syn(parser, [], by_par, par_less))
    _status_update('synonyms')
    gen_synonym_fk(engine, parser, tax_id2_tax_pk, gen_anc_first_then_syn(parser, synonyms, {}, []))
    _status_update('nothing more to do')


_col_verb_rank_norm = {None: 'infraspecies',
                       'f.': 'forma',
                       'subsp.': 'subspecies',
                       'var.': 'varietas',
                       }


def _normalize_rank(parser, rec):
    r_idx, vr_idx = parser.taxonRank_idx, parser.verbatimTaxonRank_idx
    rs, vrs = rec[r_idx], rec[vr_idx]
    if rs == 'infraspecies':
        rs = _col_verb_rank_norm[vrs]
        rec[r_idx] = rs
    else:
        assert vrs is None
    return rs


def assign_traversals(parser, by_par, par_less):
    trav_int = 0
    enter_idx = parser.entry_trav_int_idx
    exit_idx = parser.exit_trav_int_idx
    id_idx = parser.taxonID_idx
    rank_sort_number_idx = parser.rank_sort_num_idx
    internal_rec_list = []
    for rec, par_rec in gen_preorder_and_par_for_each_subtree(parser, by_par, par_less):
        rank_str = _normalize_rank(parser, rec)
        sn = RANK_TO_SORTING_NUMBER[rank_str]
        if par_rec:
            if par_rec[rank_sort_number_idx] <= sn:
                raise ValueError('Rank tension... Record:\n{}\nparent:\n{}'.format(rec, par_rec))
        rec[rank_sort_number_idx] = sn
        rec[enter_idx] = trav_int
        trav_int += 1
        rec_id = rec[id_idx]
        c_list = by_par.get(rec_id)
        if c_list:
            internal_rec_list.append(rec)
        else:
            rec[exit_idx] = rec[enter_idx]
    for rec in reversed(internal_rec_list):
        rec_id = rec[id_idx]
        c_list = by_par.get(rec_id)
        assert c_list
        rec[exit_idx] = max([i[exit_idx] for i in c_list])


def gen_dataset_fk(engine, parser, record_iter):
    _generic_gen_key(engine, record_iter=record_iter,
                     val_index=parser.datasetID_idx,
                     table_class=Dataset, val_db_name='name',
                     ini_map={'Species 2000': 0})


def gen_modified_fk(engine, parser, record_iter):
    _generic_gen_key(engine, record_iter=record_iter,
                     val_index=parser.modified_idx,
                     table_class=Modified, val_db_name='modified')


def gen_according_to_fk(engine, parser, record_iter):
    _generic_gen_key(engine, record_iter=record_iter,
                     val_index=parser.nameAccordingTo_idx,
                     table_class=AccordingTo, val_db_name='according_to')


def gen_description_fk(engine, parser, record_iter):
    _generic_gen_key(engine, record_iter=record_iter,
                     val_index=parser.description_idx,
                     table_class=Description, val_db_name='description')


def gen_taxon_concept_fk(engine, parser, record_iter):
    _generic_gen_key(engine, record_iter=record_iter,
                     val_index=parser.taxonConceptID_idx,
                     table_class=TaxonConcept, val_db_name='concept_id')


def gen_sci_name_fk(engine, parser, record_iter):
    _generic_gen_key(engine, record_iter=record_iter,
                     val_index=parser.scientificNameID_idx,
                     table_class=SciNameUID, val_db_name='sci_name_uid')


def gen_authority_fk(engine, parser, record_iter):
    _generic_gen_key(engine, record_iter=record_iter,
                     val_index=parser.scientificNameAuthorship_idx,
                     table_class=Authority, val_db_name='label')


def _gen_tsi_ini_kwarg(value):
    # These indices have to be kept in sync with the value tuples created in
    #   append_tax_sup_inf_fk
    return {'modified_fk': value[0],
            'according_to_fk': value[1],
            'description_fk': value[2],
            'taxon_concept_fk': value[3],
            'col_link_fk': value[4],
            'dataset_fk': value[5],
            'is_extinct': value[11],
            'uuid_fk': value[6],
            'sci_name_uid_fk': value[7],
            'name_status_fk': value[8],
            'syn_name_status_fk': value[9],
            }


def append_tax_sup_inf_fk(engine, parser, record_iter):
    # Taxon Supp Info gets added as separate integer to each record...
    val_to_key = {}
    next_key = 0
    tsi_idx = parser.taxon_supp_inf_idx
    for n, rec in enumerate(record_iter):
        if rec[parser.acceptedNameUsageID_idx] is None:
            ns_fk = rec[parser.taxonomicStatus_idx]
            sns_fk = None
        else:
            sns_fk = rec[parser.taxonomicStatus_idx]
            ns_fk = None
        is_extinct = rec[parser.isExtinct_idx]
        value = (rec[parser.modified_idx],  # 0
                 rec[parser.nameAccordingTo_idx],  # 1
                 rec[parser.description_idx],  # 2
                 rec[parser.taxonConceptID_idx],  # 3
                 rec[parser.references_idx],  # 4
                 rec[parser.datasetID_idx],  # 5
                 rec[parser.identifier_idx],  # 6
                 rec[parser.scientificNameID_idx],  # 7
                 ns_fk,  # 8
                 sns_fk,  # 9
                 0 if is_extinct is None else 1,  # 10
                 is_extinct,  # 11
                 )
        key = val_to_key.setdefault(value, next_key)
        if key == next_key:
            next_key += 1
        rec[tsi_idx] = key
    _add_orm_objs(engine, val_to_key, TaxonSuppInfo, kwarg_fn=_gen_tsi_ini_kwarg)


def _gen_taxon_name(value):
    ntax_words = len(value) - 1
    assert ntax_words > 0
    assert ntax_words < 5
    k = {'word_1_fk': value[0],
         'word_2_fk': None if ntax_words < 2 else value[1],
         'word_3_fk': None if ntax_words < 3 else value[2],
         'word_4_fk': None if ntax_words < 4 else value[3],
         'authority_fk': value[-1],
         }
    return k


def _gen_rank_fields(value):
    return {'rank': value[0], 'sorting_number': value[1]}


def gen_rank_fk(engine, parser, record_iter):
    rank_tuple_to_key = {}
    next_key = 0
    for n, rec in enumerate(record_iter):
        rank = rec[parser.taxonRank_idx]
        rank_num = rec[parser.rank_sort_num_idx]
        tup = (rank, rank_num)
        key = rank_tuple_to_key.setdefault(tup, next_key)
        if key == next_key:
            next_key += 1
        rec[parser.taxonRank_idx] = key
    _add_orm_objs(engine, rank_tuple_to_key, Rank, kwarg_fn=_gen_rank_fields)


def gen_name_fk(engine, parser, record_iter):
    word_tuple_to_key = {}
    word_to_key = {}
    next_word_key = 0
    next_tuple_key = 0
    for n, rec in enumerate(record_iter):
        name_parts = parse_name_parts(parser, rec, append_auth=False, include_subgenus=False)[0]
        nw = len(name_parts)
        assert nw > 0
        assert nw <= 4
        word_key_list = []
        for word in name_parts:
            key = word_to_key.setdefault(word, next_word_key)
            if key == next_word_key:
                next_word_key += 1
            word_key_list.append(key)
        word_key_list.append(rec[parser.scientificNameAuthorship_idx])
        work_key_tuple = tuple(word_key_list)
        key = word_tuple_to_key.setdefault(work_key_tuple, next_tuple_key)
        if key == next_tuple_key:
            next_tuple_key += 1
        rec[parser.scientificName_idx] = key
    _add_orm_objs(engine, word_to_key, NameFragment, val_db_name='word')
    _add_orm_objs(engine, word_tuple_to_key, TaxonName, kwarg_fn=_gen_taxon_name)


def gen_name_status_fk(engine, parser, record_iter):
    syn_to_key = {}
    next_syn_key = 0
    valid_to_key = {}
    next_valid_key = 0
    for rec in record_iter:
        v = rec[parser.acceptedNameUsageID_idx]
        name_stat = rec[parser.taxonomicStatus_idx]
        is_valid = v is None
        if is_valid:
            key = valid_to_key.setdefault(name_stat, next_valid_key)
            if next_valid_key == key:
                next_valid_key += 1
        else:
            key = syn_to_key.setdefault(name_stat, next_syn_key)
            if next_syn_key == key:
                next_syn_key += 1
        rec[parser.taxonomicStatus_idx] = key
    _add_orm_objs(engine, syn_to_key, SynonymStatus, val_db_name='label')
    _add_orm_objs(engine, valid_to_key, ValidNameStatus, val_db_name='label')


def _cache_uuid(uuid, uuid_to_key, next_uuid_key):
    uuid_key = uuid_to_key.setdefault(uuid, next_uuid_key)
    if uuid_key == next_uuid_key:
        next_uuid_key += 1
    return next_uuid_key, uuid_key


def _gen_col_link_kwargs(taxon_tup):
    # sync with indices in gen_taxon_fk
    return {'prim_uuid': taxon_tup[1],
            'syn_uuid': taxon_tup[2]
            }


def gen_uuid_and_col_link_fk(engine, parser, record_iter):
    pref = '^http://www.catalogueoflife.org/annual-checklist/2015/details/species/id/'
    uuid_pat = '([a-zA-Z0-9]+)'
    syn_pat = re.compile(r'{p}{u}/synonym/{u}$'.format(p=pref, u=uuid_pat))
    primary_pat = re.compile(r'{p}{u}$'.format(p=pref, u=uuid_pat))
    col_link_tup_to_key = {}
    col_next_key = 0
    uuid_to_key = {}
    next_uuid_key = 0
    for rec in record_iter:
        uuid = rec[parser.identifier_idx]
        next_uuid_key, uuid_key = _cache_uuid(uuid, uuid_to_key, next_uuid_key)
        rec[parser.identifier_idx] = uuid_key
        col_link = rec[parser.references_idx]
        if col_link:
            m = syn_pat.match(col_link)
            if m:
                f = m.group(1).strip()
                next_uuid_key, ffk = _cache_uuid(f, uuid_to_key, next_uuid_key)
                s = m.group(2).strip()
                next_uuid_key, sfk = _cache_uuid(s, uuid_to_key, next_uuid_key)
                col_value = (2, ffk, sfk)
            else:
                m = primary_pat.match(col_link)
                if not m:
                    raise ValueError('reference:\n{}\n did not match any pattern'.format(col_link))
                f = m.group(1).strip()
                next_uuid_key, ffk = _cache_uuid(f, uuid_to_key, next_uuid_key)
                col_value = (1, ffk, None)
        else:
            col_value = (0, None, None)
        col_key = col_link_tup_to_key.setdefault(col_value, col_next_key)
        if col_key == col_next_key:
            col_next_key += 1
        rec[parser.references_idx] = col_key
    _add_orm_objs(engine, uuid_to_key, TaxUUID, val_db_name='uuid')
    _add_orm_objs(engine, col_link_tup_to_key, CoLLink, kwarg_fn=_gen_col_link_kwargs)


def _gen_taxon_kwargs(taxon_tup):
    # sync with indices in gen_taxon_fk
    return {'name_fk': taxon_tup[0],
            'rank_fk': taxon_tup[1],
            'tax_sup_info_fk': taxon_tup[2],
            }


def _gen_edge_kwargs(edge_tup):
    # sync with indices in gen_taxon_fk
    return {'parent_fk': edge_tup[0],
            'child_fk': edge_tup[1],
            }


def gen_taxon_fk(engine, parser, record_iter):
    taxon_tup_to_key = {}
    next_key = 0
    tax_id_to_pk = {}
    edges = {}
    next_edge_key = 0
    for rec in record_iter:
        assert rec[parser.acceptedNameUsageID_idx] is None
        tup = (rec[parser.scientificName_idx],  # name
               rec[parser.taxonRank_idx],  # rank
               rec[parser.taxon_supp_inf_idx],  # tax_supp_info
               )
        key = taxon_tup_to_key.setdefault(tup, next_key)
        if next_key == key:
            next_key += 1
        rec_id = rec[parser.taxonID_idx]
        tax_id_to_pk[rec_id] = key
        par_id = rec[parser.parentNameUsageID_idx]
        if par_id is not None:
            assert par_id in tax_id_to_pk
            par_key = tax_id_to_pk[par_id]
            etup = (par_key, key)
            ekey = edges.setdefault(etup, next_edge_key)
            if ekey == next_edge_key:
                next_edge_key += 1
        del rec[:]
    _add_orm_objs(engine, taxon_tup_to_key, Taxon, kwarg_fn=_gen_taxon_kwargs)
    _add_orm_objs(engine, edges, Edge, kwarg_fn=_gen_edge_kwargs)
    return tax_id_to_pk


def _gen_syn_kwargs(taxon_tup):
    # sync with indices in gen_taxon_fk
    return {'name_fk': taxon_tup[0],
            'valid_taxon_fk': taxon_tup[1],
            'tax_sup_info_fk': taxon_tup[2],
            }


def gen_synonym_fk(engine, parser, tax_id_to_pk, record_iter):
    taxon_tup_to_key = {}
    next_key = 0
    for rec in record_iter:
        acc_tax_id = rec[parser.acceptedNameUsageID_idx]
        valid_fk = tax_id_to_pk[acc_tax_id]
        tup = (rec[parser.scientificName_idx],  # name
               valid_fk,  # valid_taxon_fk
               rec[parser.taxon_supp_inf_idx],  # tax_supp_info
               )
        key = taxon_tup_to_key.setdefault(tup, next_key)
        if next_key == key:
            next_key += 1
        del rec[:]
    _add_orm_objs(engine, taxon_tup_to_key, Synonym, kwarg_fn=_gen_syn_kwargs)


def parse_name_parts(p, rec, append_auth=False, include_subgenus=True):
    rank_str = rec[p.taxonRank_idx]
    assert rank_str
    sci_name = rec[p.scientificName_idx]
    if is_higher_taxon(rank_str):
        name_words = [sci_name]
    else:
        name_words = [rec[p.genericName_idx], rec[p.specificEpithet_idx]]
        if is_infraspecific(rank_str):
            # if not _has_printed_ssp:
            #     _has_printed_ssp = True
            #     print('a ssp: {}'.format(rec))
            # print(n, 'below species', rec)
            interjection = rec[p.verbatimTaxonRank_idx]
            if interjection:
                name_words.append(interjection.strip())
            name_words.append(rec[p.infraspecificEpithet_idx])
        # elif not _has_printed_sp:
        #     _has_printed_sp = True
        #     print('a sp: {}'.format(rec))

    subgenus = rec[p.subgenus_idx]
    ns = list(name_words)
    if len(ns) > 1 and subgenus and include_subgenus:
        ns = [ns[0]] + ['({})'.format(subgenus.strip())] + ns[1:]
    if append_auth:
        auth = rec[p.scientificNameAuthorship_idx]
        if auth:
            ns.append(auth)
    return ns, sci_name


timestamp = None
curr_operation = None


def format_td(td):
    sec_int = int(td.seconds)
    hours = sec_int // 3600
    hours = hours % 24
    minutes = (sec_int % 3600) // 60
    seconds = sec_int % 60
    return '{:2}h:{:2}m:{:2}s'.format(hours, minutes, seconds)


def _status_update(next_operation, msg=None):
    global curr_operation, timestamp
    if timestamp is not None:
        assert curr_operation is not None
        td = datetime.now() - timestamp
        if msg:
            _LOG.debug('in progress message {} (has taken {} since start)'.format(msg, format_td(td)))
            return
        else:
            _LOG.debug('In-memory processing of {} took {}'.format(curr_operation, format_td(td)))
    _LOG.debug('Starting in-memory processing of {}'.format(next_operation))
    timestamp = datetime.now()
    curr_operation = next_operation


def gen_preorder_and_par_for_each_subtree(parser, by_parent, parentless):
    n = 0
    m = '  yielding record #{:7} of in-memory operation'
    id_idx = parser.taxonID_idx
    for rec in parentless:
        next_tup = (rec, None)
        deferred_subtrees = []
        more_to_go_on_this_subtree = True
        while more_to_go_on_this_subtree:
            if next_tup is None:
                while deferred_subtrees:
                    deferred_sibs = deferred_subtrees[-1]
                    if deferred_sibs:
                        next_tup = deferred_sibs.pop(0)
                        break
                    deferred_subtrees.pop()
                if next_tup is None:
                    more_to_go_on_this_subtree = False
                    break
            n += 1
            if n % IN_MEM_STATUS_FREQ == 0:
                _status_update(None, m.format(n))
            yield next_tup
            rec_id = next_tup[0][id_idx]
            child_list = by_parent.get(rec_id)
            if child_list:
                first_child = child_list[0]
                delayed_sibs = [(c, rec) for c in child_list[1:]]
                deferred_subtrees.append(delayed_sibs)
                next_tup = (first_child, rec)
            else:
                next_tup = None


def gen_anc_first_for_each_subtree(parser, by_parent, parentless):
    n = 0
    m = '  yielding record #{:7} of in-memory operation'
    for rec in parentless:
        rec_id = rec[parser.taxonID_idx]
        to_check_for_children = [rec_id]
        n += 1
        yield rec
        while to_check_for_children:
            par_id = to_check_for_children.pop(0)
            child_list = by_parent.get(par_id)
            if child_list:
                for c in child_list:
                    rec_id = c[parser.taxonID_idx]
                    if rec_id in by_parent:
                        to_check_for_children.append(rec_id)
                    n += 1
                    if n % IN_MEM_STATUS_FREQ == 0:
                        _status_update(None, m.format(n))
                    yield c


def gen_anc_first_then_syn(parser, synonyms, by_parent, parentless):
    to_check_for_children = []
    n = 0
    m = '  yielding synon. #{:7} of in-memory operation'
    for rec in gen_anc_first_for_each_subtree(parser, by_parent, parentless):
        n += 1
        yield rec
    for r in synonyms:
        n += 1
        if n % IN_MEM_STATUS_FREQ == 0:
            _status_update(None, m.format(n))
        yield r


def _generic_gen_key(engine, val_index, record_iter, table_class, val_db_name, ini_map=None):
    _LOG.debug('replacing value for {}.{} in memory...'.format(table_class.__name__, val_db_name))
    v2k = replace_entry_with_key(record_iter, val_index=val_index, initial_mapping=ini_map)
    return _add_orm_objs(engine, v2k, table_class, val_db_name)


def _add_orm_objs(engine, v2k, table_class, val_db_name=None, kwarg_fn=None):
    vlist = list(v2k.keys())
    nv = len(vlist)
    _LOG.debug('sorting the {} value of {}.{}'.format(nv, table_class.__name__, val_db_name))
    vlist.sort(key=lambda v: (0, v) if v is None else (1, v))
    session = sessionmaker(bind=engine)()
    nr = len(v2k)
    if nr > 1000:
        if nr > 100000:
            commit_freq = 50000
        else:
            commit_freq = nr // 2
    else:
        commit_freq = max(1, nr)

    status_freq = commit_freq // 2
    n = 0
    try:
        for n, val in enumerate(vlist):
            key = v2k[val]
            if kwarg_fn is None:
                kwarg = {val_db_name: val}
            else:
                kwarg = kwarg_fn(val)
            session.add(table_class(pk=key, **kwarg))
            num = 1 + n
            if num % commit_freq == 0:
                _status_update(None, 'committing {} changes {}.{}  to the db'.format(commit_freq,
                                                                                     table_class.__name__,
                                                                                     val_db_name))
                session.commit()
            elif num % status_freq == 0:
                _status_update(None, '    adding {}/{} "{}" -> {} ...'.format(n + 1, nv, val, key))

    finally:
        m = 'committing the remaining {} of the {}.{} changes to the db'
        _status_update(None, m.format((n + 1) % commit_freq, table_class.__name__, val_db_name))
        session.commit()


IN_MEM_STATUS_FREQ = 100000


def replace_entry_with_key(record_iter, val_index, initial_mapping=None):
    if not initial_mapping:
        val_to_key = {}
        next_key = 0
    else:
        val_to_key = initial_mapping
        next_key = 1 + max(initial_mapping.values())
    for n, record in enumerate(record_iter):
        value = record[val_index]
        key = val_to_key.setdefault(value, next_key)
        if key == next_key:
            next_key += 1
        record[val_index] = key
    return val_to_key
