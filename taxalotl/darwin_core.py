#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Much of following code is from by JAR and from:
#   reference-taxonomy/feed/gbif/project_2016.py
# and
#   reference-taxonomy/feed/gbif/process_gbif_taxonomy.py
from __future__ import print_function

import io
import os
import re
# noinspection PyPep8Naming
import xml.etree.ElementTree as ET

from peyotl import (assure_dir_exists,
                    get_logger)

from taxalotl.ott_schema import InterimTaxonomyData
from taxalotl.resource_wrapper import TaxonomyWrapper
from .util import OutFile

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
_genus_re = u"[A-ZÖ{x}-]+".format(x=_lower)
_epithets_re = u"(|{e}|{e}{e}|{e}{e}{e})".format(e=_epithet)
_canon_re = u"{}{}".format(_genus_re, _epithets_re)
_auth_re = u" +(d'|von |van |de |dem |der |da |del |di |le |f\\. |[{}(])(..|\\.).*".format(_upper)
_with_subg = u"{g} \\({g}\\){e}".format(g=_genus_re, e=_epithets_re)
_trimmer = re.compile(u"({})({})".format(_canon_re, _auth_re))
_subg_trimmer = re.compile(u"({})({})".format(_with_subg, _auth_re))
_year_re = re.compile(u".*, [12][0-9][0-9][0-9?]\\)?")
_has_digit = re.compile(u".*[0-9].*")


def canonical_name(name):
    if u' phage ' in name or name.endswith(' phage'):
        return name
    if u' virus ' in name or name.endswith(' virus'):
        return name
    m = _subg_trimmer.match(name)
    if m is not None:
        canon = m.group(1)
        cs = canon.split('(')
        nc = None
        if len(cs) == 2:
            g = cs[0].strip()
            el = cs[1].split(')')
            if len(el) == 2:
                e = el[1].strip()
                nc = '{} {}'.format(g, e)
        if nc is None:
            raise RuntimeError('Could not parse name with subgenus:\n"{}"'.format(name))
        canon = nc
    else:
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


def write_gbif_projection_file(source, destination, fields2index):
    i = 0
    sci_ind = fields2index['scientificName']
    tax_id = fields2index['id']
    pnu_ind = fields2index['parentNameUsageID']
    anu_ind = fields2index['acceptedNameUsageID']
    tr_ind = fields2index['taxonRank']
    ts_ind = fields2index['taxonomicStatus']
    nat_ind = fields2index['nameAccordingTo']
    with io.open(source, 'rU', encoding='utf-8') as infile:
        with OutFile(destination) as outfile:
            for line in infile:
                row = line.split('\t')
                scientific = row[sci_ind]
                tax_id_str = row[tax_id]
                canenc = canonical_name(scientific)
                row_el = [tax_id_str,
                          row[pnu_ind],  # parentNameUsageID
                          row[anu_ind],  # acceptedNameUsageID
                          canenc,  # canonicalName
                          row[tr_ind],  # taxonRank
                          row[ts_ind],  # taxonomicStatus
                          row[nat_ind],  # nameAccordingTo / datasetID
                          ]
                row_el = [x.strip() for x in row_el]
                row_str = u"\t".join(row_el)
                outfile.write(row_str)
                outfile.write("\n")
                if i % 500000 == 0:
                    _LOG.info(u"{} {} => {}".format(i, scientific, canenc))
                i += 1


def read_gbif_projection(proj_filepath, itd, field_to_index, do_gbif_checks):
    col_taxon_id = field_to_index['id']
    col_par_name_usage_id = field_to_index['parentNameUsageID']
    col_accepted_name_usage_id = field_to_index['acceptedNameUsageID']
    col_canonical_name = field_to_index['canonicalName']
    col_taxon_rank = field_to_index['taxonRank']
    col_taxonomic_status = field_to_index['taxonomicStatus']
    col_name_according_to = field_to_index['nameAccordingTo']
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
    to_ignore = set()
    count = 0
    n_syn = 0
    with io.open(proj_filepath, 'rU', encoding='utf-8') as inp:
        for line_num, row in enumerate(inp):
            fields = row.split('\t')
            # acceptedNameUsageID
            syn_target_id_string = fields[col_accepted_name_usage_id].strip()
            is_synonym = False
            if syn_target_id_string:
                is_synonym = True
            try:
                taxon_id = int(fields[col_taxon_id])
            except:
                if line_num == 0:
                    continue
                raise
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
            if do_gbif_checks and (tstatus != 'accepted' and taxon_id not in not_doubtful):
                m = "Unexpected non accepted: id={} name=\"{}\" tstatus={} source={}"
                raise RuntimeError(m.format(taxon_id, name, tstatus, source))
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
    x = list(orphaned)
    x.sort()
    itd.details_log["orphans_pruned"] = x
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
    manifest_fp = os.path.join(source, 'meta.xml')
    manifest_root = ET.parse(manifest_fp).getroot()
    core_paths = []
    field2index = {}
    for el in manifest_root.findall('{http://rs.tdwg.org/dwc/text/}core'):
        for sub in el:
            if sub.tag.endswith('}id'):
                field2index['id'] = int(sub.attrib['index'])
            elif sub.tag.endswith('}field'):
                nns = os.path.split(sub.attrib['term'])[-1]
                field2index[nns] = int(sub.attrib['index'])
        for f in el.findall('{http://rs.tdwg.org/dwc/text/}files'):
            for loc in f.findall('{http://rs.tdwg.org/dwc/text/}location'):
                core_paths.append(loc.text.strip())
    if len(core_paths) != 1:
        raise ValueError('Did not find a single core path in DwC file ("{}") found: {}'.format(manifest_fp, core_paths))
    taxon_fn = core_paths[0]
    proj_out = os.path.join(destination, 'projection.tsv')
    if not os.path.exists(proj_out):
        proj_in = os.path.join(source, taxon_fn)
        write_gbif_projection_file(proj_in, proj_out, field2index)
    homemade = {'id': 0,
                'parentNameUsageID': 1,
                'acceptedNameUsageID': 2,
                'canonicalName': 3,
                'taxonRank': 4,
                'taxonomicStatus': 5,
                'nameAccordingTo': 6,
                }

    itd = InterimTaxonomyData()
    to_remove, to_ignore, paleos = read_gbif_projection(proj_out, itd, homemade,
                                                        do_gbif_checks=isinstance(res_wrapper, GBIFWrapper))
    add_fake_root(itd)
    remove_if_tips(itd, to_remove)
    o_to_ignore = find_orphaned(itd)
    to_ignore.update(o_to_ignore)
    prune_ignored(itd, to_ignore)
    _LOG.info('writing {} paleodb ids'.format(len(paleos)))
    with OutFile(os.path.join(destination, 'paleo.tsv')) as paleofile:
        for taxon_id in paleos:
            paleofile.write('{}\n'.format(taxon_id))
    res_wrapper.post_process_interim_tax_data(itd)
    itd.write_to_dir(destination)


_BOLD_NAME = re.compile(r"BOLD[:]([A-Z0-9]+)")
class GBIFWrapper(TaxonomyWrapper):
    schema = {"http://rs.tdwg.org/dwc/"}

    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    def normalize(self):
        normalize_darwin_core_taxonomy(self.unpacked_filepath, self.normalized_filedir, self)

    def node_should_be_semanticized(self, node):
        if _BOLD_NAME.match(node.name):
            return False
        return True