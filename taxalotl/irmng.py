#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR
from __future__ import print_function

import csv
import os
import re

from peyotl import (get_logger)

from taxalotl.ott_schema import InterimTaxonomyData
from taxalotl.resource_wrapper import TaxonomyWrapper

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
        m = "Expecting 1 file to IRMNG_DWC_SP_PROFILE.*.csv found: {}".format(poss_p)
        raise RuntimeError(m)
    if len(poss_i) == 2:
        if poss_i[0] == prof_file:
            return os.path.join(source, poss_i[1]), os.path.join(source, prof_file)
        if poss_i[1] == prof_file:
            return os.path.join(source, poss_i[1]), os.path.join(source, prof_file)
    raise RuntimeError("Expecting 2 files to match IRMNG_DWC.*.csv found: {}".format(poss_i))


def read_irmng_file(irmng_file_name):
    # 0 "TAXONID","SCIENTIFICNAME","SCIENTIFICNAMEAUTHORSHIP","GENUS",
    # 4 "SPECIFICEPITHET","FAMILY","TAXONRANK","TAXONOMICSTATUS",
    # 8 "NOMENCLATURALSTATUS","NAMEACCORDINGTO","ORIGINALNAMEUSAGEID",
    # 11 "NAMEPUBLISHEDIN","ACCEPTEDNAMEUSAGEID","PARENTNAMEUSAGE",
    # 14 "PARENTNAMEUSAGEID","TAXONREMARKS","MODIFIED","NOMENCLATURALCODE"
    itd = InterimTaxonomyData()

    rows = 0
    to_par = itd.to_par
    to_children = itd.to_children
    to_rank = itd.to_rank
    synonyms = itd.synonyms
    itd.extra_blob = {}
    to_tsta_nstat_keep = itd.extra_blob
    itd.syn_id_to_valid = {}
    syn_id_to_valid = itd.syn_id_to_valid
    with open(irmng_file_name, 'rU', encoding='utf-8') as csvfile:
        csvreader = csv.reader(csvfile)
        header = next(csvreader)
        if header[5] != 'FAMILY':
            m = 'IRMNG csv failed header check: header[5] == {} != not "FAMILY"'.format(header[5])
            raise RuntimeError(m)
        for raw_row in csvreader:
            # noinspection PyCompatibility
            row = [i for i in raw_row]
            taxon_id = int(row[0])
            long_name = row[1]
            auth = row[2]
            rank = row[6]
            tstatus = row[7]  # TAXONOMICSTATUS
            nstatus = row[8]  # NOMENCLATURALSTATUS
            try:
                syn_target_id = int(row[12]) if row[12] else None
            except:
                _LOG.warn('Dropping unparseable line {}'.format(row))
                continue
            parent = row[-4]
            diff_target = syn_target_id is not None and syn_target_id != taxon_id
            synonymp = tstatus == 'synonym' or diff_target
            # Calculate taxon name
            genus = row[3]
            if rank == 'species':
                epithet = row[4]
                name = u'{} {}'.format(genus, epithet)
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
                    synstat = nstatus if nstatus else tstatus
                    itd.register_synonym(syn_target_id, name, synstat, syn_id=taxon_id)
                    assert taxon_id != syn_target_id
                    syn_id_to_valid[taxon_id] = syn_target_id
                else:
                    _LOG.info(u"Dropping synonym without target: {} '{}'".format(taxon_id, name))
                continue
            # Kludge to get rid of redundancies e.g. Megastoma
            if tstatus == '':
                aa_found = False
                for value in row:
                    if 'awaiting allocation' in value:
                        aa_found = True
                        break
                if aa_found:
                    _LOG.info(u"Dropping awaiting allocation taxon: {} '{}'".format(taxon_id, name))
                    continue
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
            to_tsta_nstat_keep[taxon_id] = [tstatus, nstatus, False]
            rows += 1
            if rows % 250000 == 0:
                _LOG.info("{} rows {} {}".format(rows, taxon_id, name))
    _LOG.info("Processed: {} taxa, {} synonyms".format(len(to_par), len(synonyms)))
    return itd


def fix_irmng(itd):
    # Get rid of all synonym of a synonym
    syn_id_to_valid = itd.syn_id_to_valid
    to_check = set(syn_id_to_valid.keys())
    loser_synonyms = set()
    syn_check_round = 0
    taboo = set()
    while len(to_check) > len(taboo):
        to_remap = {}
        _LOG.info('syn_check_round = {}'.format(syn_check_round))
        syn_check_round += 1
        for syn_id in to_check:
            if syn_id in taboo:
                continue
            valid_id = syn_id_to_valid[syn_id]
            if valid_id in syn_id_to_valid:
                new_valid = syn_id_to_valid[valid_id]
                assert new_valid != valid_id
                to_remap[syn_id] = new_valid
        for k, v in to_remap.items():
            old_v = syn_id_to_valid[k]
            if v == k:
                _LOG.info("Synonymy ring: Arbitrarily leaving {} -> {} mapping.".format(v, old_v))
                taboo.add(k)
            else:
                itd.fix_synonym(v, old_v, k)
                syn_id_to_valid[k] = v
        to_check = set(to_remap.keys())
        loser_synonyms.update(to_check)
    ril = list(loser_synonyms)
    ril.sort()
    itd.details_log["indirect synonyms"] = ril
    _LOG.info("Indirect synonyms: {}".format(len(ril)))

    # Short-circuit taxon parents that are synonyms
    to_par = itd.to_par
    del_syn_par = set()
    par_fixes = {}
    id_to_children = itd.to_children
    for taxon_id, par_id in to_par.items():
        valid_par = syn_id_to_valid.get(par_id)
        if valid_par is not None:
            clist = id_to_children.get(par_id)
            del_syn_par.add(par_id)
            if clist is not None:
                for c_id in clist:
                    par_fixes[c_id] = valid_par
                del id_to_children[par_id]
    for taxon_id, par_id in par_fixes.items():
        to_par[taxon_id] = par_id
        id_to_children.setdefault(par_id, []).append(taxon_id)
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
    to_tsta_nstat_keep = itd.extra_blob
    keep_count = 0
    actual_grandfatherd = set()
    missing_pars = set()
    roots = itd.root_nodes
    for taxon_id in to_par.keys():
        tsta_nst_keep = to_tsta_nstat_keep[taxon_id]
        if tsta_nst_keep[2]:
            continue  # already seen
        if (taxon_id in grandfathered
                or (tsta_nst_keep[0] in taxon_statuses_to_keep
                    and tsta_nst_keep[1] in nomenclatural_statuses_to_keep)):
            scan_id = taxon_id
            while not tsta_nst_keep[2]:
                if scan_id in grandfathered:
                    _LOG.info(u'Grandfathering {}'.format(scan_id))
                    actual_grandfatherd.add(scan_id)
                tsta_nst_keep[2] = True
                keep_count += 1
                if scan_id not in to_par:
                    missing_pars.add(scan_id)
                    _LOG.info(u"Missing parents for {}".format(scan_id))
                    break
                psi = scan_id
                scan_id = to_par[scan_id]
                if not scan_id:
                    _LOG.info(u"No parent for {}".format(psi))
                    roots.add(psi)
                    break
                valid_id = scan_id
                if valid_id in syn_id_to_valid:
                    m = u'anc in syn_id_to_valid: syn_id_to_valid[{}] = {} anc of {}'
                    raise ValueError(m.format(valid_id, syn_id_to_valid[valid_id], taxon_id))
                tsta_nst_keep = to_tsta_nstat_keep.get(valid_id)
                if tsta_nst_keep is None:
                    m = u'anc ({}) of taxon {} not in to_tsta_nstat_keep'
                    _LOG.info(m.format(valid_id, taxon_id))
                    roots.add(psi)
                    break

    ids_reg = to_par.keys()
    for irmng_id in ids_reg:
        if not to_tsta_nstat_keep[irmng_id][2]:
            par_id = to_par[irmng_id]
            if par_id:
                pc = id_to_children.get(par_id)
                if pc and to_tsta_nstat_keep[irmng_id][2]:
                    pc.remove(irmng_id)
            del to_par[irmng_id]
            clist = id_to_children.get(irmng_id)
            if clist:
                del id_to_children[irmng_id]
                for child in clist:
                    if to_tsta_nstat_keep[child][2]:
                        roots.add(child)
                        to_par[child] = None

    _LOG.info("Keeping {} taxa".format(keep_count))
    _LOG.info("{} missing parents".format(len(missing_pars)))


def read_extinct_info(profile_file_name, itd):
    not_extinct = frozenset(['1531',  # Sarcopterygii
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
                             ])
    to_par = itd.to_par
    d = {}
    with open(profile_file_name, 'rU', encoding='utf-8') as csvfile:
        csvreader = csv.reader(csvfile)
        header = next(csvreader)
        if header[1] != 'ISEXTINCT':
            raise ValueError('ISEXTINCT in header row but found "{}"'.format(header[1]))
        for row in csvreader:
            taxonid = int(row[0])
            if taxonid not in to_par:
                continue
            is_extinct = (row[1] == 'TRUE')
            if taxonid in not_extinct:
                if not is_extinct:
                    _LOG.info('protected IRMNG ID {} not extinct:'.format(taxonid))
                else:
                    _LOG.info('Fixing extinctness of IRMNG ID {}:'.format(taxonid))
                    is_extinct = False
            if is_extinct:
                d[taxonid] = True
    itd.extinct_known = d


# noinspection PyUnusedLocal
def normalize_irmng(source, destination, res_wrapper):
    i_file, prof_file = _find_irmng_input_files(source)
    itd = read_irmng_file(i_file)
    fix_irmng(itd)
    read_extinct_info(prof_file, itd)
    res_wrapper.post_process_interim_tax_data(itd)
    itd.write_to_dir(destination)


class IRMNGWrapper(TaxonomyWrapper):
    schema = {"irmng dwc"}

    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    def normalize(self):
        normalize_irmng(self.unpacked_filepath, self.normalized_filedir, self)


"""
    
# Command like arguments: something like
#      feed/irmng/in/IRMNG_DWC.csv
#      feed/irmng/in/IRMNG_DWC_SP_PROFILE.csv
#      tax/irmng/taxonomy.tsv
#      tax/irmng/synonyms.tsv

import csv, sys






taxa = {}
synonyms = {}
roots = []




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

    with OutFile(taxonomy_file_name) as taxfile:
        print
        'Writing %s' % taxonomy_file_name
        taxfile.write('%s\t|\t%s\t|\t%s\t|\t%s\t|\t%s\t|\t\n' % (
        'uid', 'parent_uid', 'name', 'rank', 'flags'))
        taxfile.write('%s\t|\t%s\t|\t%s\t|\t%s\t|\t%s\t|\t\n' % ('0', '', 'life', 'no rank', ''))
        for taxon in taxa.itervalues():
            if taxon.keep:
                write_taxon(taxon, taxfile)

    with OutFile(synonyms_file_name) as synfile:
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
