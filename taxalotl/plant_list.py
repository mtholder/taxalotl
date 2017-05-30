#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR
from __future__ import print_function
from unidecode import unidecode
import csv
import codecs
import os
import time

from bs4 import BeautifulSoup as Soup
from peyotl import (assure_dir_exists, get_logger, download_large_file)

from taxalotl.resource_wrapper import TaxonomyWrapper

_LOG = get_logger(__name__)
DOMAIN = "http://www.theplantlist.org"
THROTTLE_BREAK = 10

_num_downloads_this_session = 0


def download_csv_for_family(fam_dir, fam_html_fp, url_pref):
    global _num_downloads_this_session
    fam_html_content = codecs.open(fam_html_fp, 'rU', encoding='utf-8').read()
    soup = Soup(fam_html_content, 'html.parser')
    csva = soup.find_all("a", attrs={"type": "text/csv"})
    if len(csva) != 1:
        raise RuntimeError(u"Not just 1 CSV type links in {} : {}".format(fam_html_fp, csva))
    csv_link = csva[0]
    csv_rel_url = csv_link['href']
    template = u'{}{}' if (url_pref.endswith('/') or csv_rel_url.startswith('/')) else u'{}/{}'
    fam_url = template.format(url_pref, csv_rel_url)
    fam_dest = os.path.join(fam_dir, csv_rel_url)
    if not os.path.exists(fam_dest):
        _LOG.debug(u"Starting download from url = {} to {}".format(fam_url, fam_dest))
        download_large_file(fam_url, fam_dest)
        _num_downloads_this_session += 1
        _LOG.debug(u"Download completed to .".format(fam_url, fam_dest))


def scrape_families_from_higher_group(out_dir, top_file):
    global _num_downloads_this_session
    dirname = os.path.split(top_file)[1] + '_families'
    fam_dir = os.path.join(out_dir, dirname)
    assure_dir_exists(fam_dir)
    top_content = codecs.open(top_file, 'rU', encoding='utf-8').read()
    soup = Soup(top_content, 'html.parser')
    nametree_list = soup.select("#nametree > li")
    _LOG.debug("will write to {}".format(dirname))

    for list_item in nametree_list:
        if _num_downloads_this_session != 0:
            m = "Sleeping for {} seconds to be polite to the server..."
            _LOG.debug(m.format(THROTTLE_BREAK))
            time.sleep(THROTTLE_BREAK)

        fam_link = list_item.select('a')
        assert len(fam_link) == 1
        fam_link = fam_link[0]
        fam_rel_url = fam_link['href']
        fam_name = fam_link.string.strip()
        fam_dest = os.path.join(fam_dir, fam_name + '.html')
        template = u'{}{}' if fam_rel_url.startswith('/') else u'{}/{}'
        fam_url = template.format(DOMAIN, fam_rel_url)
        if not os.path.exists(fam_dest):
            _LOG.debug(u"Starting download from url = {} to {}".format(fam_url, fam_dest))
            download_large_file(fam_url, fam_dest)
            _num_downloads_this_session += 1
            _LOG.debug(u"Download completed to .".format(fam_url, fam_dest))
        download_csv_for_family(fam_dir, fam_dest, fam_url)


def _gen_line(ls):
    try:
        return '{}\n'.format('\t|\t'.join(ls))
    except:
        conv = [unidecode(i) for i in ls]
        m = u'Could not serialize to ASCII: "{}"   converted to "{}"'
        old = u'", "'.join(ls)
        new = '", "'.join(conv)
        _LOG.warn(m.format(old, new))
        # import sys; sys.exit('{}\n{}'.format(repr(conv), repr(ls)))
        return '{}\n'.format('\t|\t'.join(conv))


pl_rank_to_ott_rank = {'f.': "form",
                       'var.': "variety",
                       'subsp.': "subspecies",
                       }
AGI = 'auto-generated-pl-id'


def normalize_plantlist_file(inp_fp, out_dir, family, maj_group_id):
    _LOG.info(u'{} to {}'.format(inp_fp, out_dir))
    fam_name = unidecode(family)
    id_to_line = {fam_name: [fam_name, maj_group_id, fam_name, 'family', AGI]}
    legit_ids = {fam_name, }
    illegit_ids = set()
    name_to_id = {}
    with codecs.open(inp_fp, 'rb') as csvfile:
        csvreader = csv.reader(csvfile)
        header = csvreader.next()
        _LOG.info(u'header = {}'.format(header))
        for n, raw_row in enumerate(csvreader):
            # noinspection PyCompatibility
            row = [unicode(i, 'utf-8') for i in raw_row]
            taxon_id = row[0]
            fam = row[2]
            if fam != family:
                raise RuntimeError("Unexpected family in taxon {} of {}: {}".format(n, family, row))
            genus = row[4]
            assert genus
            is_hybrid = bool(row[5])
            flags = 'hybrid' if is_hybrid else ''
            sp_epithet = row[6]
            infr_rank = row[7]
            infr_epi = row[8]
            par_id = None
            if infr_rank:
                rank = pl_rank_to_ott_rank[infr_rank]
                assert infr_epi
                name = ' '.join([genus, sp_epithet, infr_epi])
            else:
                if infr_epi:
                    rank = 'infraspecificname'
                    name = ' '.join([genus, sp_epithet, infr_epi])
                elif sp_epithet:
                    rank = 'species'
                    name = ' '.join([genus, sp_epithet])
                else:
                    rank = 'genus'
                    name = genus
                    par_id = fam_name
            tax_stat = row[10]
            id_to_line[taxon_id] = [taxon_id, par_id, name, rank, flags]
            if tax_stat.lower() == 'accepted':
                if name in name_to_id:
                    m = 'Name "{}" repeated in {}. IDs {} and {}. Ignoring the second...'
                    _LOG.warn(m.format(name, family, name_to_id[name], taxon_id))
                    continue
                if rank == 'species' or rank == 'genus':
                    name_to_id[name] = taxon_id
                legit_ids.add(taxon_id)
            else:
                illegit_ids.add(taxon_id)
            _LOG.info(u'taxon_id={} "{}" "{}" "{}" rank={} tax_stat={}'.format(taxon_id,
                                                                               genus,
                                                                               sp_epithet,
                                                                               infr_epi,
                                                                               rank,
                                                                               tax_stat))
    # uid	|	parent_uid	|	name	|	rank	|	flags	|
    legit_gen, legit_sp, legit_infr = [], [], []
    for vid in legit_ids:
        line_el = id_to_line[vid]
        rank = line_el[3]
        if rank in ['genus', 'family']:
            if rank != 'family':
                legit_gen.append(vid)
            par_id = line_el[1]
        elif rank == 'species':
            name = line_el[2]
            gen_name = name.split(' ')[0]
            par_id = name_to_id.get(gen_name)
            if par_id is None:
                gen_gen_id = gen_name
                assert gen_gen_id not in id_to_line
                id_to_line[gen_gen_id] = [gen_gen_id, fam_name, gen_name, 'genus', AGI]
                name_to_id[gen_name] = gen_gen_id
                legit_gen.append(gen_gen_id)
                _LOG.info("autogenerating genus record for {}".format(gen_name))
                par_id = gen_gen_id
            legit_sp.append(vid)
        else:
            name = line_el[2]
            sp_name = ' '.join(name.split(' ')[:2])
            par_id = name_to_id.get(sp_name)
            if par_id is None:
                gen_sp_id = sp_name
                assert gen_sp_id not in id_to_line
                id_to_line[gen_sp_id] = [gen_sp_id, sp_name.split()[0], sp_name, 'species', AGI]
                name_to_id[sp_name] = gen_sp_id
                _LOG.info("autogenerating species record for {}".format(sp_name))
                par_id = sp_name
                legit_sp.append(gen_sp_id)
            legit_infr.append(vid)
        line_el[1] = par_id
    id_order = legit_gen + legit_sp + legit_infr
    j = '\t|\t'
    taxon_fp = os.path.join(out_dir, 'taxonomy.tsv')
    assure_dir_exists(out_dir)
    with codecs.open(taxon_fp, 'w', encoding='utf-8') as outp:
        outp.write('{}\n'.format(j.join(['uid', 'parent_uid', 'name', 'rank', 'flags'])))
        outp.write('{}\n'.format(j.join(id_to_line[fam_name])))
        for i in id_order:
            outp.write(_gen_line(id_to_line[i]))
    not_accepted_fp = os.path.join(out_dir, 'not-accepted.tsv')
    with codecs.open(not_accepted_fp, 'w', encoding='utf-8') as outp:
        outp.write('{}\n'.format(j.join(['uid', 'name', 'rank', 'flags'])))
        for i in illegit_ids:
            line_el = id_to_line[i]
            tout = [line_el[0]] + line_el[2:]
            outp.write(_gen_line(tout))


class PlantListWrapper(TaxonomyWrapper):
    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    def download(self):
        dd = self.unpacked_filepath
        assure_dir_exists(dd)
        _LOG.info("uf = {}".format(dd))
        top_files = []
        for u in self.url_list:
            pref, suff = os.path.split(u)
            if not suff:
                pref, suff = os.path.split(pref)
            _LOG.info("p = {} s = {}".format(pref, suff))
            assert suff
            dfp = os.path.join(dd, suff)
            top_files.append(dfp)
            if not os.path.exists(dfp):
                _LOG.debug("Starting download from {} to {}".format(u, dfp))
                download_large_file(u, dfp)
                _LOG.debug("Download from {} to {} completed.".format(u, dfp))
        for dfp in top_files:
            scrape_families_from_higher_group(dd, dfp)
        open(self.download_filepath, 'w')

    def normalize(self):
        dd = self.unpacked_filepath
        subdirs = ['{}_families'.format(i) for i in 'AGPB']
        dn = self.normalized_filepath
        assure_dir_exists(dn)
        for s in subdirs:
            nsd = os.path.join(dn, s)
            assure_dir_exists(nsd)
            in_sd = os.path.join(dd, s)
            csvf = [i for i in os.listdir(in_sd) if i.endswith('.csv')]
            csvf.sort()
            for i in csvf:
                inp_fp = os.path.join(in_sd, i)
                stem = i[:-4]
                out_dir = os.path.join(nsd, stem)
                normalize_plantlist_file(inp_fp, out_dir, stem, maj_group_id=s[0])
            _LOG.info('csvf = {}'.format(stem))

        _LOG.info("dd = {} dn = {}".format(dd, subdirs))

    @property
    def download_filepath(self):
        uf = self.unpacked_filepath
        if uf is None:
            return None
        return os.path.join(uf, "download_complete.txt")
