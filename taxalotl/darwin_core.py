#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR and from:
#   reference-taxonomy/feed/gbif/project_2016.py
# and
#   reference-taxonomy/feed/gbif/process_gbif_taxonomy.py
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
                row_el = [row[1], # taxonID
                          row[3], # parentNameUsageID
                          row[4], # acceptedNameUsageID
                          canenc, # canonicalName
                          row[7], # taxonRank
                          row[10], # taxonomicStatus
                          row[2], # nameAccordingTo / datasetID
                         ]
                row_el = [i.strip() for i in row_el]
                row_str = u"\t".join(row_el)
                outfile.write(row_str)
                outfile.write("\n")
                if i % 500000 == 0:
                    _LOG.info("{} {} => {}".format(i, scientific.encode('utf-8'), canenc))
                i += 1

def read_gbif_projection(proj_filepath):
    col_taxonID = 0
    col_parentNameUsageID = 1
    col_acceptedNameUsageID = 2
    col_canonicalName = 3
    col_taxonRank = 4
    col_taxonomicStatus = 5
    col_nameAccordingTo = 6
    not_doubtful = {
        8407745: "Hierococcyx"
    }
    to_ignore = []  # stack
    to_ignore.append(0)  # kingdom incertae sedis
    flushed_because_source = set()
    to_remove = set()
    synonyms = {}
    to_par = {}  # key is the child id and the value is the parent
    to_children = {}  # key is the parent and value is the list of children
    to_rank = {}  # key is the node id and the value is the rank
    to_name = {}
    paleos = set()
    ranks_to_ignore = frozenset(["form", "variety", "subspecies", "infraspecificname"])
    to_ignore = set()
    root_nodes = set()
    count = 0
    n_syn = 0
    with codecs.open(proj_filepath, 'r', encoding='utf-8') as inp:
        for row in inp:
            fields = row.split('\t')
            # acceptedNameUsageID
            syn_target_id_string = fields[col_acceptedNameUsageID].strip()
            is_synonym = False
            if syn_target_id_string:
                is_synonym = True
            taxon_id = int(fields[col_taxonID])
            name = fields[col_canonicalName].strip()
            assert name
            source = fields[col_nameAccordingTo].strip()
            tstatus = fields[col_taxonomicStatus].strip()  # taxonomicStatus
            # Filter out IRMNG and IPNI tips,
            # See http://www.gbif.org/dataset/d9a4eedb-e985-4456-ad46-3df8472e00e8
            if (("IRMNG Homonym" in source) or
                    ("Interim Register of Marine" in source) or
                    ("International Plant Names Index" in source) or
                    (source == "d9a4eedb-e985-4456-ad46-3df8472e00e8")):
                flushed_because_source.add(taxon_id)
                if is_synonym:
                    continue
                else:
                    to_remove.add(taxon_id)
            elif is_synonym:
                synon_of = int(syn_target_id_string)
                synonyms.setdefault(synon_of, []).append((name, tstatus))
                n_syn += 1
                continue
            elif ("Paleobiology Database" in source) or (
                source == "c33ce2f2-c3cc-43a5-a380-fe4526d63650"):
                paleos.add(taxon_id)
            if tstatus == 'synonym' or (tstatus == 'doubtful' and taxon_id not in not_doubtful):
                to_remove.add(taxon_id)
                continue
            assert tstatus == 'accepted'
            rank = fields[col_taxonRank].strip()
            if rank in ranks_to_ignore:
                to_ignore.add(taxon_id)

            parent_id_string = fields[col_parentNameUsageID].strip()

            # Past all the filters, time to store
            to_name[taxon_id] = name
            to_rank[taxon_id] = rank

            if parent_id_string:
                par_id = int(parent_id_string)
                to_par[id] = par_id
                to_children.setdefault(par_id, []).append(taxon_id)
            else:
                assert rank == 'kingdom'
                to_par[taxon_id] = None
                root_nodes.add(taxon_id)

            count += 1
            if count % 100000 == 0:
                _LOG.info("lines={} #syn={} #roots={}".format(count, len(synonyms), len(root_nodes)))

def normalize_darwin_core_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    proj_out = os.path.join(destination, 'projection.tsv')
    if not os.path.exists(proj_out):
        proj_in = os.path.join(source, 'taxon.txt')
        write_gbif_projection_file(proj_in, proj_out)
    x = read_gbif_projection(proj_out)
    outfile = open(os.path.join(outdir, "taxonomy.tsv"), "w")
    outfilesy = open(os.path.join(outdir, "synonyms.tsv"), "w")

'''
def process_gbif(inpath, outdir):
    
    infile = open(inpath, "r")

    infile_taxon_count = 0
    infile_synonym_count = 0
    count = 0
    bad_id = 0
    no_parent = 0
    parent = {}  # key is taxon id, value is the parent
    children = {}  # key is taxon id, value is list of children (ids)
    nm_storage = {}  # key is taxon id, value is the name
    nrank = {}  # key is taxon id, value is rank
    synnames = {}  # key is synonym id, value is name
    syntargets = {}  # key is synonym id, value is taxon id of target
    syntypes = {}  # key is synonym id, value is synonym type
    to_remove = []  # list of ids
    paleos = []  # ids that come from paleodb
    flushed_because_source = 0

    print('%s taxa, %s synonyms\n' % (infile_taxon_count, infile_synonym_count))

    print('%s bad id; %s no parent id; %s synonyms; %s bad source' %
          (bad_id, no_parent, len(synnames), flushed_because_source))

    # Parent/child homonyms now get fixed by smasher

    # Flush terminal taxa from IRMNG and IPNI (OTT picks up IRMNG separately)
    count = 0
    for id in to_remove:
        if (not id in children):  # and id in nrank and nrank[id] != "species":
            if id in nm_storage:
                del nm_storage[id]
                # should remove from children[parent[id]] too
            count += 1
    print
    "tips removed (IRMNG and IPNI):", count

    # Put parentless taxa into the ignore list.
    # This isn't really needed any more; smasher can cope with multiple roots.
    count = 0
    for id in nm_storage:
        if id in parent and parent[id] not in nm_storage:
            count += 1
            if parent[id] != 0:
                to_ignore.append(id)
                if count % 1000 == 0:
                    print
                    "example orphan ", id, nm_storage[id]
    print
    "orphans to be pruned: ", count

    # Now delete the taxa-to-be-ignored and all of their descendants.
    if len(to_ignore) > 0:
        print
        'pruning %s taxa' % len(to_ignore)
        seen = {}
        stack = to_ignore
        while len(stack) != 0:
            curid = stack.pop()
            if curid in seen:
                continue
            seen[curid] = True
            if curid in children:
                for id in children[curid]:
                    stack.append(id)
        for id in seen:
            if id in nm_storage:
                del nm_storage[id]

    """
    output the id parentid name rank
    """
    print
    "writing %s taxa" % len(nm_storage)
    outfile.write("uid\t|\tparent_uid\t|\tname\t|\trank\t|\t\n")

    count = 0
    for id in nm_storage:
        parent_id = ""
        if id == incertae_sedis_kingdom:
            print
            "kingdom incertae sedis should have been deleted by now"
        elif id in parent:
            parent_id = str(parent[id])
        elif nrank[id] == 'kingdom':
            parent_id = "0"
        outfile.write("%s\t|\t%s\t|\t%s\t|\t%s\t|\t\n" %
                      (id, parent_id, nm_storage[id], nrank[id]))
        count += 1
        if count % 100000 == 0:
            print
            count
    outfile.write("0\t|\t\t|\tlife\t|\t\t|\t\n")
    outfile.close()

    print
    "writing %s synonyms" % len(synnames)
    outfilesy.write('uid\t|\tname\t|\ttype\t|\t\n')
    for id in synnames:
        target = syntargets[id]  # taxon id of target (int)
        if target in nm_storage:
            outfilesy.write('%s\t|\t%s\t|\t%s\t|\t\n' %
                            (target, synnames[id], syntypes[id]))
    outfilesy.close()

    print
    'writing %s paleodb ids' % len(paleos)
    paleofile = open(os.path.join(outdir, 'paleo.tsv'), 'w')
    for id in paleos:
        paleofile.write(('%s\n' % id))
    paleofile.close()



'''