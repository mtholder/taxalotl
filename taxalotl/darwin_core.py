#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR and from:
#   reference-taxonomy/feed/gbif/project_2016.py

from __future__ import print_function
import re
import codecs
import os
import time

from peyotl import (assure_dir_exists,
                    get_logger,
                    write_as_json)

_LOG = get_logger(__name__)



def normalize_darwin_core(source, destination, res_wrapper):
    pass



# Cases to deal with:
#  Foo bar
#  Foo bar Putnam
#  Foo bar Putnam, 1972
#  Foo bar Putnam, 4723     no authority
#  Foo bar Putnam 1972      no authority (in GBIF)
#  Enterobacteria phage PA-2
#  Ajuga pyramidalis L.
_lower = u"a-záåäàãçëéèïíøöóü'×?"
_upper = u"A-ZÄÁÅÁÁÇČÐĎĐÉÉÎİŁŘŠŚŞȘÔØÖÔÓÜÚŽ"
_epithet = u" +[%s0-9.-]+" % _lower

# Matches a canonical name
_canon_re = u"[A-ZÖ{l}-]+(|{e}|{e}{e}|{e}{e}{e})".format(l=_lower, e=_epithet)
_auth_re = u" +(d'|von |van |de |dem |der |da |del |di |le |f\\. |[{}(])(..|\\.).*".format(_upper)
_trimmer = re.compile(u"({})({})".format(_canon_re, _auth_re))
_year_re = re.compile(u".*, [12][0-9][0-9][0-9?]\\)?")
_has_digit = re.compile(u".*[0-9].*")

def canonical_name(name):
    if u' phage ' in name or name.endswith(' phage'):
        return name
    if u' virus ' in name or name.endswith(' virus'):
        return name
    m = _trimmer.match(name)
    if m is None:
        return name
    canon = m.group(1)
    # group 1 = canonical name
    # group 2 = epithet(s)
    # group 3 = authority
    # group 4 = capital letter or prefix
    if _has_digit.match(name):
        haz = _year_re.match(name) != None
    else:
        haz = True
    return canon if haz else name


def write_gbif_projection_file(source, destination):
    i = 0
    with codecs.open(source, 'r', encoding='utf-8') as infile:
        with codecs.open(destination, 'w', encoding='utf-8') as outfile:
            for line in infile:
                row = line.split('\t')
                scientific = row[6]
                canenc = canonical_name(scientific)
                row_el = (row[1], # taxonID
                          row[3], # parentNameUsageID
                          row[4], # acceptedNameUsageID
                          canenc, # canonicalName
                          row[7], # taxonRank
                          row[10], # taxonomicStatus
                          row[2], # nameAccordingTo / datasetID
                         )
                row_str = u"\t".join(row_el)
                outfile.write(row_str)
                outfile.write("\n")
                if i % 500000 == 0:
                    _LOG.info("{} {} => {}".format(i, scientific.encode('utf-8'), canenc))
                i += 1

def normalize_darwin_core_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    proj_out = os.path.join(destination, 'projection.tsv')
    if not os.path.exists(proj_out):
        proj_in = os.path.join(source, 'taxon.txt')
        write_gbif_projection_file(proj_in, proj_out)

