#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from peyotl import (add_or_append_to_dict, assure_dir_exists,
                    get_logger,
                    write_as_json)
import tempfile
import codecs
import os

_LOG = get_logger(__name__)


def write_ott_taxonomy_tsv(out_fp,
                           root_nodes,
                           id_to_par,
                           id_to_children,
                           id_to_rank,
                           id_to_name,
                           has_syn_dict,
                           details_log):
    """If has_syn_dict is provided, then a list of the IDs that occur in that dict
    is returned in the order that the IDs were written to taxonomy file. This 
    allows for the synonyms.tsv file to be written in a similar order, which makes browsing it
    easier.
    """
    if has_syn_dict is None:
        has_syn_dict = {}
    syn_id_order = []
    num_tips_written = 0
    num_internals_written = 0
    rn = list(root_nodes)
    rn.sort()
    with codecs.open(out_fp, 'w', encoding='utf-8') as out:
        out.write("uid\t|\tparent_uid\t|\tname\t|\trank\t|\t\n")
        # need to print id, parent id, and name
        for root_id in rn:
            stack = [root_id]
            while stack:
                curr_id = stack.pop()
                try:
                    if curr_id in has_syn_dict:
                        syn_id_order.append(curr_id)
                    name = id_to_name[curr_id]
                    par_id = id_to_par.get(curr_id)
                    if par_id is None:
                        spar_id = ''
                    else:
                        spar_id = str(par_id)
                    rank = id_to_rank.get(curr_id, '')
                    children = id_to_children.get(curr_id)
                    if children:
                        num_internals_written += 1
                        stack.extend(children)
                    else:
                        num_tips_written += 1
                    out.write(u'{}\n'.format('\t|\t'.join([str(curr_id), spar_id, name, rank, ''])))
                except:
                    _LOG.error("Error writing taxon_id {}".format(curr_id))
                    raise
    details_log['num_tips_written'] = num_tips_written
    details_log['num_internals_written'] = num_internals_written
    return syn_id_order


def write_ott_synonyms_tsv(out_fp,
                           id_to_name_name_type_list,
                           id_order,
                           details_log):
    num_syn_written = 0
    with codecs.open(out_fp, 'w', encoding='utf-8') as out:
        out.write("uid\t|\tname\t|\ttype\t|\t\n")
        for nd_id in id_order:
            syn_list = id_to_name_name_type_list[nd_id]
            for name, name_type in syn_list:
                num_syn_written += 1
                out.write(u'{}\n'.format('\t|\t'.join([str(nd_id), name, name_type, ''])))
    details_log['num_synonyms_written'] = num_syn_written
    details_log['num_ids_with_synonyms_written'] = len(id_order)


def write_ott_forwards(out_fp, forwarded_dict):
    with codecs.open(out_fp, 'w', encoding='utf-8') as out:
        for key, value in forwarded_dict.items():
            out.write('{}\t{}\n'.format(key, value))


def write_ncbi_details_json(fp, details_log):
    write_as_json(details_log, fp, indent=2)


class InterimTaxonomyData(object):
    def __init__(self):
        self.about = {}
        self.details_log = {}
        self.forwards = {}  # from old ID to new ID
        self.to_par = {}  # ID -> parent ID or None
        self.to_children = {}  # ID to list of children IDs
        self.to_rank = {}  # ID -> rank string
        self.root_nodes = set()  # set of IDs
        self.to_name = {}  # ID -> name
        self.name_to_ids = {}  # name to
        self.synonyms = {}
        self.repeated_names = set()

    def finalize(self):
        self.details_log['num_forwards'] = len(self.forwards)
        self.details_log['num_nodes'] = len(self.to_par)
        self.details_log['num_distinct_names'] = len(self.name_to_ids)
        self.details_log['num_ids_with_synonyms'] = len(self.synonyms)

    def register_id_and_name(self, taxon_id, name):
        self.to_name[taxon_id] = name
        if add_or_append_to_dict(self.name_to_ids, name, taxon_id):
            self.repeated_names.add(name)

    def register_synonym(self, valid_id, syn_name, name_type, syn_id=None):
        self.synonyms.setdefault(valid_id, []).append((syn_name, name_type))

    def write_ott_taxonomy_tsv(self, fp):
        return write_ott_taxonomy_tsv(fp,
                                      self.root_nodes,
                                      self.to_par,
                                      self.to_children,
                                      self.to_rank,
                                      self.to_name,
                                      self.synonyms,
                                      self.details_log)

    def write_to_dir(self, destination):
        # Write out in OTT form
        d = tempfile.mkdtemp()
        fn = ['taxonomy.tsv',
              'synonyms.tsv',
              'forwards.tsv',
              'about.json',
              'details.json']
        try:
            syn_order = self.write_ott_taxonomy_tsv(os.path.join(d, 'taxonomy.tsv'))
            write_ott_synonyms_tsv(os.path.join(d, 'synonyms.tsv'),
                                   self.synonyms,
                                   syn_order,
                                   self.details_log)
            if self.forwards:
                write_ott_forwards(os.path.join(d, 'forwards.tsv'), self.forwards)

            about_fp = os.path.join(d, 'about.json')
            write_as_json(self.about, about_fp, indent=2)
            self.finalize()
            write_ncbi_details_json(os.path.join(d, 'details.json'),
                                    self.details_log)
        except:
            for f in fn:
                tf = os.path.join(d, f)
                if os.path.exists(tf):
                    try:
                        os.remove(tf)
                    except:
                        pass
            try:
                os.rmdir(d)
            except:
                pass
            raise
        assure_dir_exists(destination)
        for f in fn:
            sfp = os.path.join(d, f)
            if os.path.exists(sfp):
                dfp = os.path.join(destination, f)
                os.rename(sfp, dfp)
        os.rmdir(d)

    def del_ids(self, id_list):
        to_name = self.to_name
        to_par = self.to_par
        to_children = self.to_children
        for taxon_id in id_list:
            if taxon_id in to_name:
                del to_name[taxon_id]
            if taxon_id in to_children:
                del to_children[taxon_id]
            pid = to_par.get(taxon_id)
            if pid:
                del to_par[taxon_id]
                pc = to_children.get(pid)
                try:
                    if pc:
                        pc.remove(taxon_id)
                except:
                    pass


