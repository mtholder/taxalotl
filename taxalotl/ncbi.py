from __future__ import print_function

import codecs
import os
import time

from peyotl import (assure_dir_exists,
                    get_logger,
                    write_as_json)

_LOG = get_logger(__name__)


###################################################################################################
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

####################################################################################################
# utility
def file_mod_time_to_isotime(fp):
    sec_since_epoch = os.path.getmtime(fp)
    return epoch_seconds_to_isotime(sec_since_epoch)


def epoch_seconds_to_isotime(sec_since_epoch):
    tuple_time = time.gmtime(sec_since_epoch)
    return time.strftime("%Y-%m-%dT%H:%M:%S", tuple_time)


def add_or_append_to_dict(d, k, v):
    """If dict `d` has key `k`, then the new value of that key will be appended
    onto a list containing the previous values, otherwise d[k] = v.
    Creates a lightweight multimap, but not safe if v can be None or a list.
    returns True if the k now maps to >1 value"""
    ov = d.get(k)
    if ov is None:
        d[k] = v
        return False
    if isinstance(ov, list):
        ov.append(v)
    else:
        d[k] = [ov, v]
    return True


####################################################################################################
def parse_ncbi_names_file(names_fp, details_log):
    """Takes a filepath to an NCBI names.dmp file.
    Returns tuple
         0 id_to_name: node_id int -> str
         1 names_to_ids: str -> int or [int, ...]
         2 synonyms node_id -> [(name, type of synonym))
    """
    count = 0
    id_to_name = {}
    synonyms = {}
    names_to_ids = {}
    repeated_names = set()
    with codecs.open(names_fp, "r", encoding='utf-8') as namesf:
        for line in namesf:
            # if you do \t|\t then you don't get the name class right because it is "\t|"
            spls = line.split("\t|")
            node_id = int(spls[0])
            name = spls[1].strip()
            homonc = spls[2].strip()  # can get if it is a series here
            nm_c = spls[3].strip()  # scientific name, synonym, etc.
            if "<series>" in homonc:
                name = name + " series"
            if "subgroup <" in homonc:  # corrects some nested homonyms
                name = homonc.replace("<", "").replace(">", "")
            # nm_c can hold
            # scientific name   - the name used in OTT as primary.
            # synonym
            # equivalent name  - usually misspelling or spelling variant
            # misspelling
            # authority  - always extends scientific name
            # type material  - bacterial strain as type for prokaryotic species ??
            # common name
            # genbank common name
            # blast name   - 247 of them - a kind of common name
            # in-part (e.g. Bacteria in-part: Monera)
            # includes (what polarity?)
            if nm_c == "scientific name":
                id_to_name[node_id] = name
                if add_or_append_to_dict(names_to_ids, name, node_id):
                    repeated_names.add(name)
            elif nm_c != "in-part":
                synonyms.setdefault(node_id, []).append((name, nm_c))
            count += 1
            if count % 100000 == 0:
                _LOG.info('{} lines of names'.format(count))
    _LOG.info("number of lines in names file: {}".format(count))
    _LOG.info("number of distinct scientific names: {}".format(len(names_to_ids)))
    _LOG.info("number of IDs with synonyms: {}".format(len(synonyms)))
    details_log['num_distinct_names'] = len(names_to_ids)
    details_log['num_ids_with_synonyms'] = len(synonyms)
    return id_to_name, names_to_ids, synonyms, repeated_names


def parse_ncbi_nodes_file(nodes_fp, details_log):
    """Takes a filepath to an NCBI nodes.dmp and returns 3 dict mapping an ID to:
        - parent ID (can be None or an int)
        - children list (only for internals)
        - rank string (if available
    """
    count = 0
    to_par = {}  # key is the child id and the value is the parent
    to_children = {}  # key is the parent and value is the list of children
    to_rank = {}  # key is the node id and the value is the rank
    root_nodes = []
    with codecs.open(nodes_fp, "r", encoding='utf-8') as nodesf:
        for line in nodesf:
            spls = line.split("\t|\t")
            ns = spls[0].strip()
            node_id = int(ns)
            ps = spls[1].strip()
            # oddly enough, the root node is ID 1 and has parent ID 1 in nodes.dmp
            if ps and ps != ns:
                par_id = int(ps)
            else:
                par_id = None
                root_nodes.append(node_id)
            rank = spls[2].strip()
            to_par[node_id] = par_id
            if rank:
                to_rank[node_id] = rank
            to_children.setdefault(par_id, []).append(node_id)
            count += 1
            if count % 100000 == 0:
                _LOG.info('{} lines of nodes'.format(count))
    _LOG.info("number of lines in nodes file: {}".format(count))
    details_log['num_nodes'] = len(to_par)
    return to_par, to_children, to_rank, root_nodes


def parse_ncbi_merged(fp, details_log):
    forwards_dict = {}
    if os.path.exists(fp):
        with codecs.open(fp, 'r', encoding='utf-8') as inp:
            for line in inp:
                rs = line.split('\t|')
                from_id, to_id = int(rs[0]), int(rs[1])
                forwards_dict[from_id] = to_id
    details_log['num_forwards'] = len(forwards_dict)
    _LOG.info('number of merges: {}'.format(len(forwards_dict)))
    return forwards_dict


# noinspection PyBroadException
def deal_with_adj_taxa_with_same_names(id_to_parent,
                                       id_to_children,
                                       id_to_rank,
                                       id_to_name,
                                       names_to_ids,
                                       synonyms,
                                       repeated_names,
                                       details_log):
    """Here we look for cases in which a taxon and its parent have the same name.
    1. If the parent is a genus (this is common for subgenera), then we just alter the
        child's name to the form 'NAME CHILDS-RANK NAME'
    2. If the parent is not a genus, then we remove the child ID from the tree.
    """
    suppressed_ids, renamed_ids = {}, set()
    for name in repeated_names:
        ids_with_this_name = names_to_ids[name]
        assert isinstance(ids_with_this_name, list) and len(ids_with_this_name) > 1
        adj_same_named_ids = []
        for i in ids_with_this_name:
            par_id = id_to_parent.get(i)
            if par_id and par_id in ids_with_this_name:
                insert_loc = len(adj_same_named_ids)
                for ind, el in enumerate(adj_same_named_ids):
                    if el[1] == par_id:
                        insert_loc = ind
                        break
                adj_same_named_ids.insert(insert_loc, (par_id, i))
        for par_id, child_id in adj_same_named_ids:
            pr = id_to_rank.get(par_id)
            if pr and pr.lower() == 'genus':
                # Change the child's name
                cr = id_to_rank.get(child_id, '')
                nn = '{} {} {}'.format(name, cr, name)
                assert nn not in names_to_ids  # could happen, but ugh...
                names_to_ids[nn] = child_id
                id_to_name[child_id] = nn
                try:
                    ids_with_this_name.remove(child_id)
                except:
                    pass
                renamed_ids.add(child_id)
            else:
                # suppress the child
                suppressed_ids[child_id] = par_id
                c_list = id_to_children.get(child_id, [])
                pc_list = id_to_children.get(par_id)
                pc_list.remove(child_id)
                pc_list.extend(c_list)
                for gc in c_list:
                    id_to_parent[gc] = par_id
                for syn_el in synonyms.get(child_id, []):
                    synonyms.setdefault(par_id, []).append(syn_el)
    details_log['ids_suppressed_because_same_name_as_par'] = suppressed_ids
    ril = list(renamed_ids)
    ril.sort()
    details_log['names_decorated_because_same_name_as_par'] = ril
    return suppressed_ids, renamed_ids


def deal_with_ncbi_env_samples_names(id_to_par, id_to_name, names_to_ids, details_log):
    ess = 'environmental samples'
    es_ids = names_to_ids.setdefault(ess, [])
    renamed_ids = set(es_ids)
    for es_id in es_ids:
        par_id = id_to_par[es_id]
        name = '{} {}'.format(id_to_name[par_id], ess)
        id_to_name[es_id] = name
        add_or_append_to_dict(names_to_ids, name, es_id)
    del names_to_ids[ess]
    ril = list(renamed_ids)
    ril.sort()
    details_log['names_decorated_because_env_samp'] = ril
    return renamed_ids


####################################################################################################
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
    with codecs.open(out_fp, 'w', encoding='utf-8') as out:
        out.write("uid\t|\tparent_uid\t|\tname\t|\trank\t|\t\n")
        # need to print id, parent id, and name
        for root_id in root_nodes:
            stack = [root_id]
            while stack:
                curr_id = stack.pop()
                if curr_id in has_syn_dict:
                    syn_id_order.append(curr_id)
                name = id_to_name[curr_id]
                par_id = id_to_par[curr_id]
                spar_id = str(par_id)
                rank = id_to_rank.get(curr_id, '')
                children = id_to_children.get(curr_id)
                if children:
                    num_internals_written += 1
                    stack.extend(children)
                else:
                    num_tips_written += 1
                out.write('{}\n'.format('\t|\t'.join([str(curr_id), spar_id, name, rank, ''])))
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
                out.write('{}\n'.format('\t|\t'.join([str(nd_id), name, name_type, ''])))
    details_log['num_synonyms_written'] = num_syn_written
    details_log['num_ids_with_synonyms_written'] = len(id_order)


def write_ott_forwards(out_fp, forwarded_dict):
    with codecs.open(out_fp, 'w', encoding='utf-8') as out:
        for key, value in forwarded_dict.items():
            out.write('{}\t{}\n'.format(key, value))


def write_ncbi_details_json(fp, details_log):
    write_as_json(details_log, fp, indent=2)


def normalize_ncbi(source, destination, url):
    nodes_fp = os.path.join(source, "nodes.dmp")
    names_fp = os.path.join(source, "names.dmp")
    merged_fp = os.path.join(source, "merged.dmp")
    about_obj = {"prefix": "ncbi",
                 "prefixDefinition": "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=",
                 "description": "NCBI Taxonomy",
                 "source": {"URL": url,
                            "date": file_mod_time_to_isotime(nodes_fp)
                            }
                 }
    details_log = {}
    forwarded = parse_ncbi_merged(merged_fp, details_log)
    id_to_par, id_to_children, id_to_rank, root_nodes = parse_ncbi_nodes_file(nodes_fp,
                                                                              details_log)
    # Make sure there is only 1 root, and that its parent is an empty string
    assert len(root_nodes) == 1
    root_id = root_nodes[0]
    assert id_to_par[root_id] is None
    id_to_par[root_id] = ""
    id_to_name, name_to_ids, synonyms, repeated_names = parse_ncbi_names_file(names_fp,
                                                                              details_log)
    # Change the root's name from root to life
    assert id_to_name[root_id] == "root"
    id_to_name[root_id] = "life"
    # Clean up adjacent taxa with the same name
    # TODO: do we really want to suppress any? Seems dangerous... MTH Apr, 2017
    deal_with_adj_taxa_with_same_names(id_to_par,
                                       id_to_children,
                                       id_to_rank,
                                       id_to_name,
                                       name_to_ids,
                                       synonyms,
                                       repeated_names,
                                       details_log)
    # Decorate the names of environmental samples
    deal_with_ncbi_env_samples_names(id_to_par, id_to_name, name_to_ids, details_log)
    # Write out in OTT form
    assure_dir_exists(destination)
    syn_order = write_ott_taxonomy_tsv(os.path.join(destination, 'taxonomy.tsv'),
                                       root_nodes,
                                       id_to_par,
                                       id_to_children,
                                       id_to_rank,
                                       id_to_name,
                                       synonyms,
                                       details_log)
    write_ott_synonyms_tsv(os.path.join(destination, 'synonyms.tsv'),
                           synonyms,
                           syn_order,
                           details_log)
    write_ncbi_details_json(os.path.join(destination, 'details.json'),
                            details_log)
    write_ott_forwards(os.path.join(destination, 'forwards.tsv'), forwarded)

    about_fp = os.path.join(destination, 'about.json')
    write_as_json(about_obj, about_fp, indent=2)

###################################################################################################
