from __future__ import print_function

import io
import os
import time
import logging
from peyutil import add_or_append_to_dict
from ..ott_schema import InterimTaxonomyData
from ..resource_wrapper import TaxonomyWrapper

_LOG = logging.getLogger(__name__)


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
#  - remove "unclassified"GENUS_RANK_TO_SORTING_NUMBER
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


####################################################################################################
def parse_ncbi_names_file(names_fp, itd):
    """Takes a filepath to an NCBI names.dmp file.
    Returns tuple
         0 id_to_name: node_id int -> str
         1 names_to_ids: str -> int or [int, ...]
         2 synonyms node_id -> [(name, type of synonym))
    """
    count = 0
    with io.open(names_fp, "r", encoding="utf-8") as namesf:
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
                itd.register_id_and_name(node_id, name)
            elif nm_c != "in-part":
                itd.register_synonym(valid_id=node_id, syn_name=name, name_type=nm_c)
            count += 1
            if count % 100000 == 0:
                _LOG.info("{} lines of names".format(count))
    _LOG.info("number of lines in names file: {}".format(count))
    _LOG.info("number of distinct scientific names: {}".format(len(itd.name_to_ids)))
    _LOG.info("number of IDs with synonyms: {}".format(len(itd.synonyms)))


def parse_ncbi_nodes_file(nodes_fp, itd):
    """Takes a filepath to an NCBI nodes.dmp and returns 3 dict mapping an ID to:
    - parent ID (can be None or an int)
    - children list (only for internals)
    - rank string (if available
    """
    count = 0
    to_par = itd.to_par  # key is the child id and the value is the parent
    to_children = itd.to_children  # key is the parent and value is the list of children
    to_rank = itd.to_rank  # key is the node id and the value is the rank
    root_nodes = itd.root_nodes
    with io.open(nodes_fp, "r", encoding="utf-8") as nodesf:
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
                root_nodes.add(node_id)
            rank = spls[2].strip()
            to_par[node_id] = par_id
            if rank:
                to_rank[node_id] = rank
            to_children.setdefault(par_id, []).append(node_id)
            count += 1
            if count % 100000 == 0:
                _LOG.info("{} lines of nodes".format(count))
    _LOG.info("number of lines in nodes file: {}".format(count))


def parse_ncbi_merged(fp, itd):
    if os.path.exists(fp):
        with io.open(fp, "r", encoding="utf-8") as inp:
            for line in inp:
                rs = line.split("\t|")
                from_id, to_id = int(rs[0]), int(rs[1])
                itd.forwards[from_id] = to_id
    _LOG.info("number of merges: {}".format(len(itd.forwards)))


# noinspection PyBroadException
def deal_with_adj_taxa_with_same_names(itd):
    """Here we look for cases in which a taxon and its parent have the same name.
    1. If the parent is a genus (this is common for subgenera), then we just alter the
        child's name to the form 'NAME CHILDS-RANK NAME'
    2. If the parent is not a genus, then we remove the child ID from the tree.
    """
    id_to_parent = itd.to_par
    id_to_children = itd.to_children
    id_to_rank = itd.to_rank
    id_to_name = itd.to_name
    names_to_ids = itd.name_to_ids
    synonyms = itd.synonyms
    repeated_names = itd.repeated_names
    suppressed_ids = {}
    renamed_ids = set()
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
            if pr and pr.lower() == "genus":
                # Change the child's name
                cr = id_to_rank.get(child_id, "")
                nn = "{} {} {}".format(name, cr, name)
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
    itd.details_log["ids_suppressed_because_same_name_as_par"] = suppressed_ids
    ril = list(renamed_ids)
    ril.sort()
    itd.details_log["names_decorated_because_same_name_as_par"] = ril


def deal_with_ncbi_env_samples_names(itd):
    id_to_par = itd.to_par
    id_to_name = itd.to_name
    names_to_ids = itd.name_to_ids
    ess = "environmental samples"
    es_ids = names_to_ids.setdefault(ess, [])
    renamed_ids = set(es_ids)
    for es_id in es_ids:
        par_id = id_to_par[es_id]
        name = "{} {}".format(id_to_name[par_id], ess)
        id_to_name[es_id] = name
        add_or_append_to_dict(names_to_ids, name, es_id)
    del names_to_ids[ess]
    ril = list(renamed_ids)
    ril.sort()
    itd.details_log["names_decorated_because_env_samp"] = ril
    return renamed_ids


####################################################################################################


def normalize_ncbi(source, destination, res_wrapper):
    url = res_wrapper.url
    itd = InterimTaxonomyData()
    nodes_fp = os.path.join(source, "nodes.dmp")
    names_fp = os.path.join(source, "names.dmp")
    merged_fp = os.path.join(source, "merged.dmp")
    itd.about = {
        "prefix": "ncbi",
        "prefixDefinition": "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=",
        "description": "NCBI Taxonomy",
        "source": {"URL": url, "date": file_mod_time_to_isotime(nodes_fp)},
    }
    parse_ncbi_merged(merged_fp, itd)
    parse_ncbi_nodes_file(nodes_fp, itd)
    # Make sure there is only 1 root, and that its parent is an empty string
    assert len(itd.root_nodes) == 1
    root_id = list(itd.root_nodes)[0]
    assert itd.to_par[root_id] is None
    itd.to_par[root_id] = ""
    parse_ncbi_names_file(names_fp, itd)
    # Change the root's name from root to life
    assert itd.to_name[root_id] == "root"
    itd.to_name[root_id] = "life"
    # Clean up adjacent taxa with the same name
    # TODO: do we really want to suppress any? Seems dangerous... MTH Apr, 2017
    deal_with_adj_taxa_with_same_names(itd)
    # Decorate the names of environmental samples
    deal_with_ncbi_env_samples_names(itd)
    res_wrapper.post_process_interim_tax_data(itd)
    itd.write_to_dir(destination)


###################################################################################################


class NCBIWrapper(TaxonomyWrapper):
    schema = {"ncbi taxonomy"}

    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    def normalize(self):
        normalize_ncbi(self.unpacked_filepath, self.normalized_filedir, self)

    def _post_process_tree(self, tree):
        self.collapse_incertae_sedis_by_name_prefix(tree, "unclassified ")

    def post_process_interim_tax_data(self, interim_tax_data):
        self.collapse_as_incertae_sedis_interim_tax_data(
            interim_tax_data, "unclassified"
        )
