#!/usr/bin/env python3

import os
import sys
from datetime import datetime
from peyotl import get_logger
from sqlalchemy import create_engine
from taxalotl.col_ingest.col_db_model import Base

from taxalotl.col_ingest.col_parser import CoLDumpParser, ParseCondenser
from taxalotl.col_ingest.key_generator import replace_values_with_keys, parse_name_parts

begin_ingest_datetime = None
NUM_SYN_INGESTED, NUM_VALID_INGESTED, FIRST_INGESTED = 0, 0, None
NUM_SYN_REC, NUM_VALID_REC = 0, 0

_LOG = get_logger(__name__)


def get_engine():
    try:
        pg_user = os.environ['PGUSER']
        pg_port = os.environ['PGPORT']
        pg_password = os.environ['PGPASSWORD']
        pg_host = os.environ['PGHOST']
    except Exception:
        sys.exit("""Sorry, still in dev mode...
    Currently you have to have PGUSER, PGPORT, PGPASSWORD, and PGHOST in your environment.
    """)
    engine_protocol = 'postgresql://{}:{}@{}:{}'.format(pg_user, pg_password, pg_host, pg_port)
    engine = create_engine('{}/col'.format(engine_protocol))
    return engine


def _validate_record(p, n, rec):
    ns, sci_name = parse_name_parts(p, rec, append_auth=True)
    expected_name = ' '.join(ns)
    if expected_name != sci_name:
        sys.exit('name unexpected: "{}" != "{}"\n{}'.format(expected_name, sci_name, rec))


def load_into_db(engine, parser, synonyms, by_parent, parentless):
    global begin_ingest_datetime, FIRST_INGESTED
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)
    replace_values_with_keys(engine, parser, synonyms, by_parent, parentless)

    """
    
    next_output_inc = 50000
    next_output = next_output_inc
    nrl = len(parentless)
    to_check_for_children = load_record_list(session, parser, parentless)
    while to_check_for_children:
        par_id = to_check_for_children.pop(0)
        child_list = by_parent.get(par_id)
        if child_list:
            nrl += len(child_list)
            to_check_for_children.extend(load_record_list(session, parser, child_list))
            if nrl > next_output:
                _LOG.debug('loaded {} records'.format(nrl))
                next_output += next_output_inc
            del by_parent[par_id]
    if True:
        _LOG.warn("Skipping synonyms!!!")
        return
    begin_ingest_datetime, FIRST_INGESTED = None, None
    _LOG.debug('to load {} synonyms...'.format(len(synonyms)))
    x = load_record_list(session, parser, synonyms)
    assert len(x) == 0
    """


def _emit_time_estimate(tot_rec, ingested, type_tag, rec_id, name):
    global begin_ingest_datetime
    pref = 'ingesting {} {:7d} of {}: {} "{}" ...'
    pref = pref.format(type_tag, ingested, tot_rec, rec_id, name)
    if begin_ingest_datetime is None:
        suffix = ''
    else:
        ts = datetime.now()
        td = ts - begin_ingest_datetime
        rtsec = td.seconds
        num_ingested = ingested - FIRST_INGESTED
        rate = num_ingested / float(rtsec)
        sec_int = int(float(tot_rec - num_ingested) / rate)
        hours = sec_int // 3600
        days = hours // 24
        hours = hours % 24
        minutes = (sec_int % 3600) // 60
        seconds = sec_int % 60
        suffix = 'so far: {:6}s@{:.2f}/sec. time left: {}d;{}h:{:2}m:{:2}s'
        suffix = suffix.format(rtsec, rate, days, hours, minutes, seconds)
    _LOG.debug('{:90}{}'.format(pref, suffix))


"""
def load_record_list(session, parser, record_list):
    global NUM_SYN_INGESTED, NUM_VALID_INGESTED, begin_ingest_datetime, FIRST_INGESTED
    STAT_OUTPUT_FREQ = 100
    if not record_list:
        return []
    # m = 'load_record_list with {} records in list. first is {}'
    # _LOG.debug(m.format(len(record_list), record_list[0][parser.scientificName_idx]))
    to_check_for_children = []
    try:
        for rec in record_list:
            if os.path.exists('STOP.txt'):
                raise RuntimeError('STOP.txt file found')
            in_db, is_valid_taxon, rec_id = record_is_already_in_db(session, rec, parser,
                                                                    to_check_for_children)
            if is_valid_taxon:
                NUM_VALID_INGESTED += 1
                if NUM_VALID_INGESTED % STAT_OUTPUT_FREQ == 0:
                    _emit_time_estimate(NUM_VALID_REC, NUM_VALID_INGESTED, 'valid taxon',
                                        rec_id, rec[parser.scientificName_idx])
            else:
                NUM_SYN_INGESTED += 1
                if NUM_SYN_INGESTED % STAT_OUTPUT_FREQ == 0:
                    _emit_time_estimate(NUM_SYN_REC, NUM_SYN_INGESTED, 'synonym',
                                        rec_id, rec[parser.scientificName_idx])
            if in_db:
                continue
            if begin_ingest_datetime is None:
                FIRST_INGESTED = NUM_VALID_INGESTED - 1
                begin_ingest_datetime = datetime.now()
            taxname_fk = _ingest_taxname_id(session, rec, parser)
            tax_sup_info_fk = _ingest_tax_sup_info(session, rec, parser)
            is_valid, stat_fk, rec_id = analyze_status(session, rec, parser)
            if is_valid:
                to_check_for_children.append(rec_id)
                _ingest_taxon(session, rec, parser, rec_id, taxname_fk, tax_sup_info_fk, stat_fk)
            else:
                _ingest_synonym(session, rec, parser, rec_id, taxname_fk, tax_sup_info_fk, stat_fk)
    finally:
        session.commit()
    return to_check_for_children


"""


def main(col_export_dir):
    global NUM_SYN_REC, NUM_VALID_REC
    if os.path.exists('STOP.txt'):
        raise RuntimeError('STOP.txt file found')

    cdp = CoLDumpParser(col_export_dir)
    p = ParseCondenser(cdp)

    engine = get_engine()
    dropped_records = []
    synonyms = []
    parentless = []
    by_parent = {}
    for n, rec in enumerate(p.gen_records()):
        try:
            _validate_record(p, n, rec)
            accepted = rec[p.acceptedNameUsageID_idx]
            if accepted is not None:
                synonyms.append(rec)
            else:
                NUM_VALID_REC += 1
                par_ind = rec[p.parentNameUsageID_idx]

                if par_ind is None:
                    parentless.append(rec)
                else:
                    by_parent.setdefault(par_ind, []).append(rec)
            if (n + 1) % 100000 == 0:
                _LOG.debug('reading {:7d}... "{}" ...'.format(n + 1, rec[p.scientificName_idx]))
        except Exception:
            _LOG.warn('Error parsing line {}:\n{}'.format(n + 1, rec))
            dropped_records.append(rec)
    NUM_SYN_REC = len(synonyms)
    load_into_db(engine, p, synonyms, by_parent, parentless)


if __name__ == '__main__':
    try:
        assert os.path.isdir(sys.argv[1])
    except Exception:
        sys.exit('Expecting the first argument to be the directory holding an CoL export.\n')
    main(sys.argv[1])
