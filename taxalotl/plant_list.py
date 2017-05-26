#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR
from __future__ import print_function
import os
import codecs
import time
from peyotl import (assure_dir_exists, get_logger, download_large_file)

from taxalotl.ott_schema import InterimTaxonomyData
from taxalotl.resource_wrapper import TaxonomyWrapper
import time
from bs4 import BeautifulSoup as Soup
_LOG = get_logger(__name__)
DOMAIN = "http://www.theplantlist.org"
THROTTLE_BREAK = 30

_num_downloads_this_session = 0

def download_csv_for_family(fam_dir, fam_html_fp, url_pref):
    global _num_downloads_this_session
    fam_html_content = codecs.open(fam_html_fp, 'rU', encoding='utf-8').read()
    soup = Soup(fam_html_content, 'html.parser')
    csva = soup.find_all("a", attrs={"type": "text/csv"})
    if len(csva) != 1:
        raise RuntimeError("Not just 1 CSV type links in {} : {}".format(fam_html_fp, csva))
    csv_link = csva[0]
    csv_rel_url = csv_link['href']
    template = '{}{}' if (url_pref.endswith('/') or csv_rel_url.startswith('/')) else '{}/{}'
    fam_url = template.format(url_pref, csv_rel_url)
    fam_dest = os.path.join(fam_dir, csv_rel_url)
    if not os.path.exists(fam_dest):
        _LOG.debug("Starting download from url = {} to {}".format(fam_url, fam_dest))
        download_large_file(fam_url, fam_dest)
        _num_downloads_this_session += 1
        _LOG.debug("Download completed to .".format(fam_url, fam_dest))

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
        if  _num_downloads_this_session != 0:
            m = "Sleeping for {} seconds to be polite to the server..."
            _LOG.debug(m.format(THROTTLE_BREAK))
            time.sleep(THROTTLE_BREAK)

        fam_link = list_item.select('a')
        assert len(fam_link) == 1
        fam_link = fam_link[0]
        fam_rel_url = fam_link['href']
        fam_name = fam_link.string.strip()
        url_pref = fam_name
        fam_dest = os.path.join(fam_dir, fam_name + '.html')
        template = '{}{}' if fam_rel_url.startswith('/') else '{}/{}'
        fam_url = template.format(DOMAIN, fam_rel_url)
        if not os.path.exists(fam_dest):
            _LOG.debug("Starting download from url = {} to {}".format(fam_url, fam_dest))
            download_large_file(fam_url, fam_dest)
            _num_downloads_this_session += 1
            _LOG.debug("Download completed to .".format(fam_url, fam_dest))
        download_csv_for_family(fam_dir, fam_dest, fam_url)



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

    @property
    def download_filepath(self):
        return os.path.join(self.unpacked_filepath, "download_complete.txt")
