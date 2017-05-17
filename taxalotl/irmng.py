#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR and from:
#   reference-taxonomy/feed/gbif/project_2016.py
# and
#   reference-taxonomy/feed/gbif/process_gbif_taxonomy.py
from __future__ import print_function

import csv
import os
import re

from peyotl import (get_logger)

from taxalotl.interim_taxonomy_struct import InterimTaxonomyData

_LOG = get_logger(__name__)


def _find_irmng_input_files(source):
    irmng_file_pat = re.compile(r"IRMNG_DWC.*\.csv")
    irmng_profile_pat = re.compile(r"IRMNG_DWC_SP_PROFILE.*\.csv")
    files = os.listdir(source)
    poss_i = [i for i in files if irmng_file_pat.match(i)]
    poss_p = [i for i in poss_i if irmng_profile_pat.match(i)]
    if len(poss_p) == 1:
        prof_file = poss_p[0]
    else:
        m = "Expecting 1 file to IRMNG_DWC_SP_PROFILE.*\.csv found: {}".format(poss_p)
        raise RuntimeError(m)
    if len(poss_i) == 2:
        if poss_i[0] == prof_file:
            return os.path.join(source, poss_i[1]), os.path.join(source, prof_file)
        if poss_i[1] == prof_file:
            return os.path.join(source, poss_i[1]), os.path.join(source, prof_file)
    raise RuntimeError("Expecting 2 files to match IRMNG_DWC.*\.csv found: {}".format(poss_i))


def read_irmng_file(irmng_file_name):
    # 0 "TAXONID","SCIENTIFICNAME","SCIENTIFICNAMEAUTHORSHIP","GENUS",
    # 4 "SPECIFICEPITHET","FAMILY","TAXONRANK","TAXONOMICSTATUS",
    # 8 "NOMENCLATURALSTATUS","NAMEACCORDINGTO","ORIGINALNAMEUSAGEID",
    # 11 "NAMEPUBLISHEDIN","ACCEPTEDNAMEUSAGEID","PARENTNAMEUSAGE",
    # 14 "PARENTNAMEUSAGEID","TAXONREMARKS","MODIFIED","NOMENCLATURALCODE"
    itd = InterimTaxonomyData()

    rows = 0
    unallocated = set()
    to_par = itd.to_par
    to_children = itd.to_children
    to_rank = itd.to_rank
    synonyms = itd.synonyms
    itd.extra_blob = {}
    to_tsta_nstat_keep_exinct = itd.extra_blob
    itd.syn_id_to_valid = {}
    syn_id_to_valid = itd.syn_id_to_valid
    with open(irmng_file_name, 'rb') as csvfile:
        csvreader = csv.reader(csvfile)
        header = csvreader.next()
        if header[5] != 'FAMILY':
            m = 'IRMNG csv failed header check: header[5] == {} != not "FAMILY"'.format(header[5])
            raise RuntimeError(m)
        for row in csvreader:
            taxon_id = int(row[0])
            long_name = row[1]
            auth = row[2]
            rank = row[6]
            tstatus = row[7]  # TAXONOMICSTATUS
            nstatus = row[8]  # NOMENCLATURALSTATUS
            syn_target_id = int(row[12]) if row[12] else None
            parent = row[-4]
            diff_target = syn_target_id is not None and syn_target_id != taxon_id
            synonymp = tstatus == 'synonym' and diff_target
            # Calculate taxon name
            genus = row[3]
            if rank == 'species':
                epithet = row[4]
                name = '{} {}'.format(genus, epithet)
            elif rank == 'genus':
                name = genus
            elif rank == 'family':
                family = row[5]
                name = family
            elif auth and long_name.endswith(auth):
                name = long_name[-len(auth) - 1:]
            else:
                name = long_name
            if synonymp:
                if diff_target:
                    itd.register_synonym(syn_target_id, name, tstatus, syn_id=taxon_id)
                    syn_id_to_valid[taxon_id] = syn_target_id
                    continue
            # Kludge to get rid of redundancies e.g. Megastoma
            if tstatus == '':
                for value in row:
                    if 'awaiting allocation' in value:
                        tstatus = 'lose'
                        unallocated.add(taxon_id)
                        break
            if parent == '':
                itd.root_nodes.add(taxon_id)
                parent = None
            else:
                parent = int(parent)

            to_par[taxon_id] = parent
            itd.register_id_and_name(taxon_id, name)
            if parent:
                to_children.setdefault(parent, []).append(taxon_id)
            to_rank[taxon_id] = rank
            to_tsta_nstat_keep_exinct[taxon_id] = [tstatus, nstatus, False, False]
            rows += 1
            if rows % 250000 == 0:
                _LOG.info("{} rows {} {}".format(rows, taxon_id, name))
    _LOG.info("Processed: {} taxa, {} synonyms".format(len(to_par), len(synonyms)))
    to_tsta_nstat_keep_exinct[101163]

    _LOG.info("Flushing: {} unallocated".format(len(unallocated)))
    return itd


def fix_irmng(itd):
    # Get rid of all synonym of a synonym
    synonyms = itd.synonyms
    syn_id_to_valid = itd.syn_id_to_valid
    loser_synonyms = set()
    for syn_id, valid_id in syn_id_to_valid.items():
        if valid_id in syn_id_to_valid:
            loser_synonyms.add(syn_id)
    ril = list(loser_synonyms)
    ril.sort()
    itd.details_log["indirect synonyms"] = ril
    _LOG.info("Indirect synonyms: {}".format(len(ril)))

    # Short-circuit taxon parents that are synonyms
    to_par = itd.to_par
    del_syn_par = set()
    par_fixes = {}
    for taxon_id, par_id in to_par.items():
        if par_id in syn_id_to_valid:
            anc_id = par_id
            while anc_id in syn_id_to_valid and anc_id not in del_syn_par:
                assert anc_id not in to_par
                del_syn_par.add(anc_id)
                anc_id = syn_id_to_valid[anc_id]
            par_fixes[taxon_id] = anc_id
    for taxon_id, par_id in par_fixes.items():
        to_par[taxon_id] = par_id
    ril = list(del_syn_par)
    ril.sort()
    itd.details_log["Parental synonyms"] = ril
    _LOG.info("Parental synonyms: {}".format(len(ril)))

    # Decide which taxa to keep
    nomenclatural_statuses_to_keep = frozenset(['',
                                                'conservandum',
                                                'nom. cons.',
                                                'Nom. cons.',
                                                'nom. cons. des.',
                                                'nom. cons. prop.',
                                                'protectum',
                                                'Nomen novum',
                                                'correctum',
                                                'orth. cons.',
                                                'legitimate',
                                                'www.nearctica.com/nomina/beetle/colteneb.htm',
                                                'Ruhberg et al., 1988',
                                                'later usage',
                                                ])
    grandfathered = frozenset([10180190, 11899420, 11704707, 10527330, 11399158, 10527966, 11444963,
                               11078615, 10522666, 10692084, 10525002, 10520170, 11444785, 11167068,
                               10531957, 11024850, 11078603, 11458858, 11081142, 11390044, 10793056,
                               10525092, 10692824, 10689467, 10543655, 10530648, 102843, 10697026,
                               10184114, 11256401, 11083597, 11102182, 11103647, 10532033, 1435408,
                               10532250, 10537012, 1178867, 10532020, 1407317, 10957072])
    taxon_statuses_to_keep = frozenset(['accepted', 'valid', ''])
    # Decide which taxa to keep
    to_tsta_nstat_keep_exinct = itd.extra_blob
    to_tsta_nstat_keep_exinct[101163]
    keep_count = 0
    actual_grandfatherd = set()
    missing_pars = set()
    for taxon_id in to_par.keys():
        tsta_nst_keep_extinct = to_tsta_nstat_keep_exinct[taxon_id]
        if tsta_nst_keep_extinct[2]:
            continue  # already seen
        if (taxon_id in grandfathered
            or (tsta_nst_keep_extinct[0] in taxon_statuses_to_keep
                and tsta_nst_keep_extinct[1] in nomenclatural_statuses_to_keep)):
            scan_id = taxon_id
            while not tsta_nst_keep_extinct[2]:
                if scan_id in grandfathered:
                    _LOG.info('Grandfathering {}', scan_id)
                    actual_grandfatherd.add(scan_id)
                tsta_nst_keep_extinct[2] = True
                keep_count += 1
                if scan_id not in to_par:
                    missing_pars.add(scan_id)
                    _LOG.info("Missing parents for {}".format(scan_id))
                    break
                psi = scan_id
                scan_id = to_par[scan_id]
                if not scan_id:
                    _LOG.info("None parent for {}".format(psi))
                    break
                valid_id = scan_id
                while valid_id in syn_id_to_valid:
                    _LOG.info('valid_id in syn_id_to_valid = {}'.format(valid_id))
                    valid_id = syn_id_to_valid[valid_id]
                tsta_nst_keep_extinct = to_tsta_nstat_keep_exinct[valid_id]

    # Now we delete the "loser synonyms"
    for syn_id in loser_synonyms:
        v = syn_id_to_valid[syn_id]
        del syn_id_to_valid[syn_id]
        if v in synonyms:
            del synonyms[v]

    _LOG.info("Keeping {} taxa".format(keep_count))
    _LOG.info("{} missing parents".format(len(missing_pars)))
    """
    # Read the file that has the extinct annotations

    with open(profile_file_name, 'rb') as csvfile:
        csvreader = csv.reader(csvfile)
        header = csvreader.next()
        if header[1] != 'ISEXTINCT':
            print >> sys.stderr, "** Expected to find ISEXTINCT in header row but didn't:", header[
                1]
        for row in csvreader:
            taxonid = row[0]
            taxon = taxa.get(taxonid)
            if taxon == None: continue
            taxon.extinctp = (row[1] == 'TRUE')
            if taxonid in not_extinct:
                if not taxon.extinctp:
                    print >> sys.stderr, 'Already not extinct: %s(%s)' % (taxonid, taxon.name)
                else:
                    print >> sys.stderr, 'Fixing extinctness of %s(%s)' % (taxonid, taxon.name)
                    taxon.extinctp = False
    """


def normalize_irmng(source, destination, res_wrapper):
    i_file, prof_file = _find_irmng_input_files(source)
    itd = read_irmng_file(i_file)
    fix_irmng(itd)
    # extinctness_report()
    # write_irmng()


"""
    
# Command like arguments: something like
#      feed/irmng/in/IRMNG_DWC.csv
#      feed/irmng/in/IRMNG_DWC_SP_PROFILE.csv
#      tax/irmng/taxonomy.tsv
#      tax/irmng/synonyms.tsv

import csv, sys


not_extinct = ['1531',  # Sarcopterygii
               '10565',  # Saurischia
               '118547',  # Aviculariidae
               '1402700',  # Trophomera
               # '11919',    # Didelphimorphia
               # '1021564',  # Cruciplacolithus
               # '1530',     # Actinopterygii

               # '1170022',  # Tipuloidea
               # '1340611',  # Retaria
               # '1124871',  # Labyrinthulomycetes [Labyrinthomorpha??]
               # '102024',   # Ophiurinidae - problem is Ophiurina
               # '1064058',  # Rhynchonelloidea genus/superfamily
               # '1114655',  # Tetrasphaera - different from GBIF
               ]






taxa = {}
synonyms = {}
roots = []


class Taxon:
    def __init__(self, id, parentid, name, rank, tstatus, nstatus):
        self.id = id
        self.parentid = parentid
        self.name = name
        self.rank = rank
        self.tstatus = tstatus
        self.nstatus = nstatus
        self.keep = False
        self.extinctp = False



def extinctness_report():
    # Report on nonextinct descended from extinct

    count = 0
    for taxon in taxa.itervalues():
        if taxon.keep and not taxon.extinctp:
            parentid = taxon.parentid
            parent = taxa.get(parentid)
            if parent != None and parent.extinctp:
                count += 1
                if taxon.rank != 'species':
                    print >> sys.stderr, ("Extant taxon %s(%s) with extinct parent %s(%s)" %
                                          (taxon.id, taxon.name, parentid, parent.name))

    print >> sys.stderr, 'Extant taxa with extinct parent:', count


# Write it out

# Returns True if one or more children also got written

def write_irmng():
    def write_taxon(taxon, taxfile):
        parentid = taxon.parentid
        if parentid == '':
            parentid = '0'
        flags = ''
        if taxon.extinctp:
            flags = 'extinct'
        taxfile.write('%s\t|\t%s\t|\t%s\t|\t%s\t|\t%s\t|\t\n' % (
        taxon.id, parentid, taxon.name, taxon.rank, flags))

    with open(taxonomy_file_name, 'w') as taxfile:
        print
        'Writing %s' % taxonomy_file_name
        taxfile.write('%s\t|\t%s\t|\t%s\t|\t%s\t|\t%s\t|\t\n' % (
        'uid', 'parent_uid', 'name', 'rank', 'flags'))
        taxfile.write('%s\t|\t%s\t|\t%s\t|\t%s\t|\t%s\t|\t\n' % ('0', '', 'life', 'no rank', ''))
        for taxon in taxa.itervalues():
            if taxon.keep:
                write_taxon(taxon, taxfile)

    with open(synonyms_file_name, 'w') as synfile:
        print
        'Writing %s' % synonyms_file_name
        synfile.write('uid\t|\tname\t|\ttype\t|\t\n')
        for syn in synonyms.itervalues():
            taxon = taxa.get(syn.parentid)
            if taxon != None and taxon.keep and not taxon.extinctp:
                status = syn.nstatus
                if status == '':
                    status = syn.tstatus
                    if status == '': status = 'synonym'
                synfile.write('%s\t|\t%s\t|\t%s\t|\t\n' % (syn.parentid, syn.name, status.lower()))



"""
