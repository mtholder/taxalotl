#!/usr/bin/env python
from __future__ import print_function

import copy
import os
import logging
import json

from peyutil import read_as_json

from ..tax_partition import (
    INP_TAXONOMY_DIRNAME,
    MISC_DIRNAME,
    get_taxon_partition,
    use_tax_partitions,
)

_LOG = logging.getLogger(__name__)
_LIFE = "Life"
####################################################################################################
# Some data (to later be refactored
_x = {
    "Archaea": {},
    "Bacteria": {},
    "Eukaryota": {
        "Archaeplastida": {
            "Glaucophyta": {},
            "Rhodophyta": {},
            "Chloroplastida": {},
            MISC_DIRNAME: {},
        },
        "Fungi": {},
        "Haptophyta": {},
        "Metazoa": {
            "Annelida": {},
            "Arthropoda": {
                "Arachnida": {},
                "Malacostraca": {},
                "Insecta": {
                    "Diptera": {},
                    "Coleoptera": {},
                    "Lepidoptera": {},
                    "Hymenoptera": {},
                    MISC_DIRNAME: {},
                },
                MISC_DIRNAME: {},
            },
            "Bryozoa": {},
            "Chordata": {},
            "Cnidaria": {},
            "Ctenophora": {},
            "Mollusca": {},
            "Nematoda": {},
            "Platyhelminthes": {},
            "Porifera": {},
            MISC_DIRNAME: {},
        },
        "SAR": {},
        MISC_DIRNAME: {},
    },
    "Viruses": {},
    MISC_DIRNAME: {},
}
BASE_PARTITIONS_DICT = {_LIFE: _x}
del _x
NAME_TO_PARTS_SUBSETS = {}
NAME_TO_PARENT_FRAGMENT = {}
NONTERMINAL_PART_NAMES = []
TERMINAL_PART_NAMES = []
PART_NAME_TO_FRAGMENT = {}


def get_inp_taxdir(parts_dir, frag, taxonomy_id):
    return os.path.join(parts_dir, frag, INP_TAXONOMY_DIRNAME, taxonomy_id)


def get_misc_inp_taxdir(parts_dir, frag, taxonomy_id):
    return os.path.join(
        parts_dir, frag, MISC_DIRNAME, INP_TAXONOMY_DIRNAME, taxonomy_id
    )


def _fill_parts_indices(d, par_frag):
    global NAME_TO_PARTS_SUBSETS, NAME_TO_PARENT_FRAGMENT, NONTERMINAL_PART_NAMES
    for k, subd in d.items():
        NAME_TO_PARTS_SUBSETS[k] = tuple(subd.keys())
        NAME_TO_PARENT_FRAGMENT[k] = par_frag
        if par_frag:
            cf = os.path.join(par_frag, k)
        else:
            cf = k
        PART_NAME_TO_FRAGMENT[k] = cf
        if subd:
            NONTERMINAL_PART_NAMES.append(k)
            _fill_parts_indices(subd, cf)
        elif k != MISC_DIRNAME:
            TERMINAL_PART_NAMES.append(k)


_fill_parts_indices(BASE_PARTITIONS_DICT, "")
PART_NAMES = list(NAME_TO_PARTS_SUBSETS.keys())
PART_NAMES.sort()
PART_NAMES = tuple(PART_NAMES)
PREORDER_PART_LIST = tuple(NONTERMINAL_PART_NAMES)
# POSTORDER_PART_LIST = tuple(reversed(PREORDER_PART_LIST))
NONTERMINAL_PART_NAMES.sort()
NONTERMINAL_PART_NAMES = tuple(NONTERMINAL_PART_NAMES)
TERMINAL_PART_NAMES.sort()
TERMINAL_PART_NAMES = tuple(TERMINAL_PART_NAMES)


def _rec_populate(d_to_fill, key_to_filled_set):
    # _LOG.info('key_to_filled_set = {}'.format(key_to_filled_set))
    if MISC_DIRNAME in d_to_fill:
        del d_to_fill[MISC_DIRNAME]
    for key, subd in d_to_fill.items():
        filled_set = key_to_filled_set.get(key)
        if subd:
            _rec_populate(subd, key_to_filled_set)
            if not filled_set:
                cu = set()
                for v in subd.keys():
                    fsv = key_to_filled_set.get(v)
                    if fsv:
                        cu.update(fsv)
                if cu:
                    key_to_filled_set[key] = cu


# Data above here, to be refactored at some point
####################################################################################################
# Code below
def iter_existing_tax_dirs(path_pref, res_id):
    suffix = os.path.join(INP_TAXONOMY_DIRNAME, res_id)
    misc_suffix = os.path.join(MISC_DIRNAME, INP_TAXONOMY_DIRNAME, res_id)
    for tup in os.walk(path_pref):
        dirname = tup[0]
        if dirname == path_pref:
            continue
        p = os.path.join(dirname, suffix)
        if os.path.exists(p):
            yield p
        p = os.path.join(dirname, misc_suffix)
        if os.path.exists(p):
            yield p


def has_any_partition_dirs(path_pref, res_id):
    assert path_pref
    for p in iter_existing_tax_dirs(path_pref, res_id):
        return True
    return False


def get_all_partition_dirs(path_pref, res_id):
    assert path_pref
    return list(iter_existing_tax_dirs(path_pref, res_id))


def find_partition_dirs_for_taxonomy(path_pref, res_id):
    return [i for i in iter_existing_tax_dirs(path_pref, res_id)]


def write_info_for_res(outstream, res, part_name_to_split):
    _LOG.debug(f"part_name_to_split = {part_name_to_split}")
    par_frag = NAME_TO_PARENT_FRAGMENT[part_name_to_split]
    _LOG.debug("par_frag = {}".format(par_frag))
    if par_frag and not res.has_been_partitioned_for_fragment(par_frag):
        par_name = os.path.split(par_frag)[-1]
        outstream.write(
            "{} does not cover or has not been partitioned into {}\n".format(
                res.id, par_name
            )
        )
        return
    part_keys = NAME_TO_PARTS_SUBSETS[part_name_to_split]
    master_map = res.get_primary_partition_map()
    mapping = [(k, master_map[k]) for k in part_keys if k in master_map]
    if not mapping:
        outstream.write("No {} mapping for {}\n".format(res.id, part_name_to_split))
        return
    fragment = (
        os.path.join(par_frag, part_name_to_split) if par_frag else part_name_to_split
    )
    if not res.has_been_partitioned_for_fragment(fragment):
        outstream.write(
            "{} does not cover or has not been partitioned into {}\n".format(
                res.id, fragment
            )
        )
        return
    outstream.write(
        "What can I say about {} at {} ? Great stuff...\n".format(res.id, fragment)
    )


def do_partition(res, strategy, part_name_to_split):
    """Partition a parent taxon into descendants and garbage (__misc__) dir

    :param res: a wrapper around the resource. Used for id, part_source_filepath,
    :param strategy: "hard-coded", "previous", or dynamic
    :param part_name_to_split must be one of the hard-coded keys in NAME_TO_PARENT_FRAGMENT
    """
    if strategy == "hard-coded":
        return do_hard_coded_partition(res, part_name_to_split)
    if strategy == "previous":
        return do_partition_from_previous(res, part_name_to_split)
    raise NotImplementedError("dynamic partitioning.")


def do_partition_from_previous(res, part_name_to_split):
    taxalotl_config = res._config
    ott = taxalotl_config.get_terminalized_res_by_id("ott", "")
    part_root_name_blob = ott.get_part_clade_names_and_blobs()
    print(res.__dict__)
    # import sys; sys.exit(json.dumps(part_root_name_blob, indent=2))
    raise NotImplementedError("previous strategy")


def do_hard_coded_partition(res, part_name_to_split):
    """Partition a parent taxon into descendants and garbage (__misc__) dir

    :param res: a wrapper around the resource. Used for id, part_source_filepath,
    :param part_name_to_split must be one of the hard-coded keys in NAME_TO_PARENT_FRAGMENT
    """
    _LOG.debug(f"part_name_to_split = {part_name_to_split}")
    par_frag = NAME_TO_PARENT_FRAGMENT[part_name_to_split]
    _LOG.debug(f"par_frag = {repr(par_frag)}")
    if par_frag and not res.has_been_partitioned_for_fragment(par_frag):
        par_name = os.path.split(par_frag)[-1]
        do_partition(res, hard_coded=True, part_name_to_split=par_name)
    part_keys = NAME_TO_PARTS_SUBSETS[part_name_to_split]
    _LOG.debug(f"part_keys = {part_keys}")
    master_map = res.get_primary_partition_map()
    _LOG.debug(f"master_map = {master_map}")
    mapping = [(k, master_map[k]) for k in part_keys if k in master_map]
    _LOG.debug(f"mapping = {mapping}")
    if not mapping:
        _LOG.info("No {} mapping for {}".format(res.id, part_name_to_split))
        return
    fragment = (
        os.path.join(par_frag, part_name_to_split) if par_frag else part_name_to_split
    )
    _LOG.debug(f"fragment = {fragment}")
    if res.has_been_partitioned_for_fragment(fragment):
        _LOG.info("Partition for fragment {} has already been done.".format(fragment))
        return
    tp = get_taxon_partition(res, fragment)
    if not par_frag:
        tp.external_input_fp = os.path.join(res.partition_source_dir, "taxonomy.tsv")
    tp.do_partition(mapping)
