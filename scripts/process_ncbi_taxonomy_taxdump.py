#!/usr/bin/env python
from __future__ import print_function
# Author: Stephen Smith

# Arguments:
#   download - T or F - whether or not to download the tar.gz file from NCBI
#   downloaddir - where to find (or put) the tar.gz and its contents
#   kill list file
#   destination dir - where taxonomy.tsv etc. are to be put

# JAR copied this file from data/ in the taxomachine repository
# to smasher/ in the opentree repository on 2013-04-25.
# Some subsequent modifications:
#  - remove "unclassified"
#  - add command line argument for directory in which to put ncbi
#  - change skipids from list to dictionary for speed
import sys
import os
from taxalotl import ResourceManager, TaxalotlConfig
from peyotl.utility import assure_dir_exists, download_large_file
from collections import Counter

"""
this processes the ncbi taxonomy tables for the synonyms and the 
names that will be included in the upload to the taxomachine
"""

"""
skipping
- X 
-environmental
-unknown
-unidentified
-endophyte
-uncultured
-scgc
-libraries
-unclassifed

if it is a series based on names 3rd column
adding series to the name

connecting these to there parents
-unclassified
-incertae sedis
"""


def download_and_unpack(url, parent_dir):
    pass


def main(raw_ncbi_dir, out_dir, download_url, download=False):
    assure_dir_exists(out_dir)
    if download:
        download_and_unpack(download_url, raw_ncbi_dir)
    print(raw_ncbi_dir, download)


if __name__ == "__main__":
    import argparse

    description = "Writes an Open Tree interim taxonomy formatted version fo NCBI's taxonomy"
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--download", action="store_true", default=False,
                   help="If present, the ncbi taxonomy will be downloaded to the \"raw\" directory")
    p.add_argument("--raw-dir", type=str,
                   help="Path to that is the parent of the raw NCBI download (default from config file)")
    p.add_argument("--out-dir", type=str,
                   help="Path to will hold the OT form of the taxonomy (default from config)")
    p.add_argument("--url", type=str, help="URL to download the raw taxonomy from.")
    p.add_argument("--resources-dir", type=str, help="the resources directory (optional)")
    p.add_argument("--config", type=str, help="the taxalotl.conf filepath (optional)")
    args = p.parse_args()
    tc = TaxalotlConfig(filepath=args.config, resources_dir=args.resources_dir)
    if args.resources_dir:
        cfg = ResourceManager(args.resources_dir)
        if not args.url:
            args.url = cfg.get('ncbi_url')
        if not args.raw_dir:
            args.raw_dir = cfg.get('ncbi_raw_dir')
        if not args.out_dir:
            args.out_dir = cfg.get('ncbi_ot_formatted_dir')
    if args.download and not args.url:
        sys.exit("If --download is used, then --url must be provided.")
    main(raw_ncbi_dir=args.raw_dir,
         out_dir=args.out_dir,
         download_url=args.url,
         download=args.download)

'''
    outfile = open(taxdir+"/taxonomy.tsv","w")
    outfilesy = open(taxdir+"/synonyms.tsv","w")

    outfile.write("uid\t|\tparent_uid\t|\tname\t|\trank\t|\t\n")

    #need to print id, parent id, and name   
    for i in lines:
        spls = lines[i].split("\t|\t")
        node_id = spls[0].strip()
        prid = parent_ids[spls[0]].strip()
        sname = spls[1].strip()

        #changed from sname to nm_storage to fix the dup name issue
        if i in final_nm_storage:
            nametowrite = final_nm_storage[i]
        else:
            nametowrite = nm_storage[i]

        # if it is the root node then we need to make its parent id blank and rename it "life"
        if nametowrite == "root":
            nametowrite = "life"
            prid = ""
        elif nametowrite == 'environmental samples':
            nametowrite = nm_storage[parent_ids[i]] + ' ' + nametowrite
            if False:
                if nametowrite not in synonyms:
                    synonyms[i] = []
                # kludge, would be better to change synonyms table representation
                synonyms[i].append('%s\t|\t%s\t|\t%s\t|\t%s\t|\t\n' % (i, 'environmental samples', '', 'synonym'))
        rankwrite = nrank[node_id]
        outfile.write(node_id+"\t|\t"+prid+"\t|\t"+nametowrite+"\t|\t"+rankwrite+"\t|\t\n")

    outfile.close()

    outfilesy.write("uid\t|\tname\t|\ttype\t|\t\n")

    for i in synonyms:
        if i in lines:
            for j in synonyms[i]:
                spls = j.split("\t|\t")
                node_id = spls[0].strip()
                sname = spls[1].strip()
                nametp = spls[3].strip()
                outfilesy.write(node_id+"\t|\t"+sname+"\t|\t"+nametp+"\t|\t\n")
    outfilesy.close()

    mergedfilename = downloaddir + "/merged.dmp"
    if os.path.isfile(mergedfilename):
        merge_count = 0
        with open(mergedfilename, 'r') as mergedfile:
            with open(taxdir + '/forwards.tsv', 'w') as forwardsfile:
                for line in mergedfile:
                    row = line.split('|')
                    from_id = row[0].strip()
                    to_id = row[1].strip()
                    forwardsfile.write("%s\t%s\n" % (from_id, to_id))
                    merge_count += 1
        print 'number of merges:', merge_count
'''
