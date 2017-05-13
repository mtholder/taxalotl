#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR and from:
#   reference-taxonomy/feed/gbif/project_2016.py
# and
#   reference-taxonomy/feed/gbif/process_gbif_taxonomy.py
from __future__ import print_function

import codecs
import os
import re

from peyotl import (assure_dir_exists,
                    get_logger)

from taxalotl.interim_taxonomy_struct import InterimTaxonomyData

_LOG = get_logger(__name__)

"""
import argparse, csv

def _read_ncbi_accession_to_ncbi_taxon_id(fp):
    att_dict = {}
    i = 0
    with codecs.open(fp, 'r', encoding='utf-8') as ato_file:
        for line in ato_file:
            ls = line.strip()
            if not ls:
                continue
            fields = line.split('\t')
            taxon_id_str = fields[1]
            if taxon_id_str != '*':
                acc = fields[0].strip()
                if len(fields) > 2 and fields[2]:
                    strain = fields[2].strip()
                else:
                    strain = None
                att_dict[acc] = (int(taxon_id_str), strain)
                i += 1
                if i % 100000 == 0:
                    _LOG.info('{} {} {}'.format(acc, taxon_id_str, strain))
    return att_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='NCBI ids and names')
    parser.add_argument('ncbi', help='taxonomy.tsv for NCBI taxonomy in open tree format')
    parser.add_argument('mappings', help='genbank accession id to NCBI taxon mapping')
    parser.add_argument('out', help='where to write the output file')
    args = parser.parse_args()
    a2t_filepath = args.mappings
    accession_to_taxon = _read_ncbi_accession_to_ncbi_taxon_id(a2t_filepath)

    print 'reading', args.ncbi
    ncbi_to_name = get_ncbi_to_name(args.ncbi)

    print 'writing', args.out
    with open(args.out, 'w') as outfile:
        i = 0
        writer = csv.writer(outfile, delimiter='\t')
        for gid in accession_to_taxon:
            (ncbi_id, strain) = accession_to_taxon[gid]
            name = ncbi_to_name.get(ncbi_id)
            row = writer.writerow((gid, ncbi_id, strain, name))
            i += 1
            if i % 500000 == 0:
                print gid, ncbi_id, strain, name
"""
def normalize_silva_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    itd = InterimTaxonomyData()
    sys.exit('hi from normalize silva')