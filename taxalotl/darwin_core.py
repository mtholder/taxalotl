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

from taxalotl.ott_schema import InterimTaxonomyData
from taxalotl.resource_wrapper import ExternalTaxonomyWrapper

_LOG = get_logger(__name__)

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
        haz = _year_re.match(name) is not None
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
                row_el = [row[1],  # taxonID
                          row[3],  # parentNameUsageID
                          row[4],  # acceptedNameUsageID
                          canenc,  # canonicalName
                          row[7],  # taxonRank
                          row[10],  # taxonomicStatus
                          row[2],  # nameAccordingTo / datasetID
                          ]
                row_el = [x.strip() for x in row_el]
                row_str = u"\t".join(row_el)
                outfile.write(row_str)
                outfile.write("\n")
                if i % 500000 == 0:
                    _LOG.info(u"{} {} => {}".format(i, scientific, canenc))
                i += 1


def read_gbif_projection(proj_filepath, itd):
    col_taxon_id = 0
    col_par_name_usage_id = 1
    col_accepted_name_usage_id = 2
    col_canonical_name = 3
    col_taxon_rank = 4
    col_taxonomic_status = 5
    col_name_according_to = 6
    not_doubtful = {
        8407745: "Hierococcyx"
    }
    flushed_because_source = set()
    to_remove = set()
    to_par = itd.to_par
    to_children = itd.to_children
    to_rank = itd.to_rank
    paleos = set()
    ranks_to_ignore = frozenset(["form", "variety", "subspecies", "infraspecificname"])
    # kingdom incertae sedis is 0
    to_ignore = {0}
    count = 0
    n_syn = 0
    with codecs.open(proj_filepath, 'r', encoding='utf-8') as inp:
        for row in inp:
            fields = row.split('\t')
            # acceptedNameUsageID
            syn_target_id_string = fields[col_accepted_name_usage_id].strip()
            is_synonym = False
            if syn_target_id_string:
                is_synonym = True
            taxon_id = int(fields[col_taxon_id])
            name = fields[col_canonical_name].strip()
            assert name
            source = fields[col_name_according_to].strip()
            tstatus = fields[col_taxonomic_status].strip()  # taxonomicStatus
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
                itd.register_synonym(synon_of, name, tstatus)
                n_syn += 1
                continue
            elif ("Paleobiology Database" in source) or (
                        source == "c33ce2f2-c3cc-43a5-a380-fe4526d63650"):
                paleos.add(taxon_id)
            if tstatus == 'synonym' or (tstatus == 'doubtful' and taxon_id not in not_doubtful):
                to_remove.add(taxon_id)
                continue
            if tstatus != 'accepted' and taxon_id not in not_doubtful:
                m = "Unexpected non accepted: {} {} {} {}".format(taxon_id,
                                                                  name,
                                                                  tstatus,
                                                                  source)
                raise RuntimeError(m)
            rank = fields[col_taxon_rank].strip()
            if rank in ranks_to_ignore:
                to_ignore.add(taxon_id)

            parent_id_string = fields[col_par_name_usage_id].strip()

            # Past all the filters, time to store
            itd.register_id_and_name(taxon_id, name)
            to_rank[taxon_id] = rank
            if parent_id_string:
                par_id = int(parent_id_string)
                to_par[taxon_id] = par_id
                to_children.setdefault(par_id, []).append(taxon_id)
            else:
                assert rank == 'kingdom'
                to_par[taxon_id] = None
                itd.root_nodes.add(taxon_id)

            count += 1
            if count % 100000 == 0:
                _LOG.info("lines={} #syn={} #roots={}".format(count,
                                                              len(itd.synonyms),
                                                              len(itd.root_nodes)))
    ril = list(flushed_because_source)
    ril.sort()
    itd.details_log["ids_suppressed_based_on_source"] = ril
    itd.details_log["num_paleodb_ids"] = len(paleos)
    ril = list(ranks_to_ignore)
    ril.sort()
    itd.details_log["ranks_ignored"] = ril
    ril = list(not_doubtful.keys())
    ril.sort()
    itd.details_log["ids_read_as_doubtful_but_retained"] = ril
    itd.details_log["removed_sources"] = ["Interim Register of Marine",
                                          "IRMNG Homonym",
                                          "International Plant Names Index",
                                          "d9a4eedb-e985-4456-ad46-3df8472e00e8"]
    return to_remove, to_ignore, paleos


def remove_if_tips(itd, to_remove):
    to_children = itd.to_children
    to_del = [i for i in to_remove if not to_children.get(i)]
    count = len(to_del)
    itd.del_ids(to_del)
    _LOG.info("tips removed (IRMNG and IPNI or status): {}".format(count))
    itd.details_log["tips_removed_because_of_source_or_status"] = count


def find_orphaned(itd):
    orphaned = set()
    to_name = itd.to_name
    id_list = to_name.keys()
    to_par = itd.to_par
    for taxon_id in id_list:
        pid = to_par.get(taxon_id)
        if pid and (pid not in to_name):
            orphaned.add(taxon_id)
    _LOG.info("orphans to be pruned: {}".format(len(orphaned)))
    l = list(orphaned)
    l.sort()
    itd.details_log["orphans_pruned"] = l
    return orphaned


def prune_ignored(itd, to_ignore):
    # Now delete the taxa-to-be-ignored and all of their descendants.
    _LOG.info('pruning {} taxa'.format(len(to_ignore)))
    seen = set()
    stack = list(to_ignore)
    stack.sort()
    itd.details_log['ignore_prune_ids'] = list(stack)
    to_children = itd.to_children
    while stack:
        curid = stack.pop()
        if curid in seen:
            continue
        seen.add(curid)
        for cid in to_children.get(curid, []):
            stack.append(cid)
    itd.del_ids(seen)


def add_fake_root(itd):
    itd.to_children[0] = list(itd.root_nodes)
    itd.to_name[0] = "life"
    itd.to_par[0] = None
    for i in itd.root_nodes:
        itd.to_par[i] = 0
    itd.root_nodes = {0}


# noinspection PyUnusedLocal
def normalize_darwin_core_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    proj_out = os.path.join(destination, 'projection.tsv')
    if not os.path.exists(proj_out):
        proj_in = os.path.join(source, 'taxon.txt')
        write_gbif_projection_file(proj_in, proj_out)
    itd = InterimTaxonomyData()
    to_remove, to_ignore, paleos = read_gbif_projection(proj_out, itd)
    add_fake_root(itd)
    remove_if_tips(itd, to_remove)
    o_to_ignore = find_orphaned(itd)
    to_ignore.update(o_to_ignore)
    prune_ignored(itd, to_ignore)
    _LOG.info('writing {} paleodb ids'.format(len(paleos)))
    with open(os.path.join(destination, 'paleo.tsv'), 'w') as paleofile:
        for taxon_id in paleos:
            paleofile.write('{}\n'.format(taxon_id))
    itd.write_to_dir(destination)


class GBIFWrapper(ExternalTaxonomyWrapper):
    def __init__(self, obj, parent=None, refs=None):
        ExternalTaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    def normalize(self):
        normalize_darwin_core_taxonomy(self.unpacked_filepath, self.normalized_filepath, self)
