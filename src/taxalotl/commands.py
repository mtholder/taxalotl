#!/usr/bin/env python
from __future__ import print_function

import os
import subprocess
import sys

from peyutil import read_as_json, write_as_json
from peyotl import (
    filter_otifacts_by_type,
    partition_otifacts_by_root_element,
    read_all_otifacts,
)
from .cmds.partitions import (
    do_partition,
    GEN_MAPPING_FILENAME,
    NAME_TO_PARTS_SUBSETS,
    PART_NAMES,
    PREORDER_PART_LIST,
    TERMINAL_PART_NAMES,
    write_info_for_res,
)
from .tax_partition import (
    INP_TAXONOMY_DIRNAME,
    MISC_DIRNAME,
    use_tax_partitions,
)


# from .cmds.analyze_update import analyze_update_to_resources
from .util import unlink, VirtCommand, OutFile
import logging

_LOG = logging.getLogger(__name__)
out_stream = sys.stdout

SEP_NAMES = "__separator_names__.json"
SEP_MAPPING = "__separator_names_to_dir__.json"
NEW_SEP_FILENAME = "__sep__.json"


# def analyze_update(taxalotl_config, id_list, level_list):
#     assert len(id_list) == 2
#     eid, lid = id_list
#     earlier = taxalotl_config.get_terminalized_res_by_id(eid)
#     later = taxalotl_config.get_terminalized_res_by_id(lid)
#     if earlier.base_id != later.base_id:
#         m = 'Can only analyze updates of the same taxonomy base: {}( base = {}), but {} (base = {})'
#         raise ValueError(m.format(eid, earlier.base_id, lid, later.base_id))
#     analyze_update_to_resources(taxalotl_config, earlier, later, level_list)


def download_resources(taxalotl_config, id_list):
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, "download")
        lic_urls, lic_tou = taxalotl_config.get_known_license_info(rw)
        if lic_urls or lic_tou:
            prompt = "The following license-related URLs were found:\n  "
            prompt += "\n  ".join(lic_urls)
            prompt += (
                "\n\nThe following license or terms of use information was found:\n\n  "
            )
            prompt += "\n\n  ".join(lic_tou)
        else:
            prompt = "No stored license info based on the last pull from the OTifacts repository."
        tag = (
            "\n\nOther license or terms of use may apply if the OTifacts repository is not "
            "up-to-date.\nEnter y to continue downloading {} : ".format(rid)
        )
        prompt += tag
        _LOG.info(repr(prompt))
        try:
            # noinspection PyCompatibility
            resp = input(prompt)
        except NameError:
            resp = input(prompt)
        if resp != "y":
            _LOG.info(
                "download of {} skipped due to lack of affirmative response.".format(
                    rid
                )
            )
        else:
            if rw.has_been_downloaded():
                m = "{} was already present at {}"
                _LOG.info(m.format(rw.id, rw.download_filepath))
            else:
                rw.download()


def _group_by_status(res, id_list):
    nd_list = []
    dnu_list = []
    unn_list = []
    n_list = []
    a_list = []
    p_list = []
    for i in id_list:
        r = res[i]
        if r.is_abstract:
            a_list.append(i)
        elif r.has_been_partitioned():
            p_list.append(i)
        elif r.has_been_normalized():
            n_list.append(i)
        elif r.has_been_unpacked():
            unn_list.append(i)
        elif r.has_been_downloaded():
            dnu_list.append(i)
        else:
            nd_list.append(i)
    return [
        ["abstract classes", a_list],
        ["not downloaded", nd_list],
        ["downloaded, but not unpacked", dnu_list],
        ["unpacked, but not normalized", unn_list],
        ["normalized", n_list],
        ["parititioned", p_list],
    ]


def get_list_of_all_resources(taxalotl_config):
    res = taxalotl_config.resources_mgr.resources
    id_list = list(res.keys())
    id_list.sort()
    if not id_list:
        m = """taxalotl does not know about any resources. This means that the resources directory
    is empty (or its contents are unparse-able).
You probably need to run the pull-otifacts command. If that does NOT solve the problem should
    report this bug and try moving that directory (so that taxalotl will create a clean one).
"""
        out_stream.write(m)
    return id_list


def status_of_resources(
    taxalotl_config, id_list, ids_only=False, by_status=False, terminal_only=False
):
    terminalize = True
    if not id_list:
        id_list = get_list_of_all_resources(taxalotl_config)
        if not id_list:
            return
        terminalize = False
    if terminal_only:
        x = []
        for i in id_list:
            ri = taxalotl_config.get_terminalized_res_by_id(i, "")
            if ri.id == i:
                x.append(i)
        id_list = x
    res = taxalotl_config.resources_mgr.resources
    if by_status:
        t_and_id_list = _group_by_status(res, id_list)
    else:
        t_and_id_list = [["", id_list]]
    # correct a wart in which "ott" is separated from its version numbers by "ott-id-list"
    tmp = []
    for tag, id_list in t_and_id_list:
        mod_id_list = list(id_list)
        if "ott" in mod_id_list and "ott-id-list" in mod_id_list:
            mod_id_list.remove("ott")
            last_ind = 0
            for index, el in enumerate(mod_id_list):
                if el.startswith("ott-id-list"):
                    last_ind = index
            mod_id_list.insert(1 + last_ind, "ott")
        tmp.append([tag, mod_id_list])
    t_and_id_list = tmp
    # End wart correction
    if ids_only:
        for tag, id_list in t_and_id_list:
            if tag:
                pref = "{}: ".format(tag)
                sep = " "
            else:
                pref, sep = "", "\n"
            if id_list:
                out_stream.write("{}{}\n".format(pref, sep.join(id_list)))
        return
    par_id_set = set()
    written = set()
    for tag, id_list in t_and_id_list:
        if tag:
            out_stream.write("{}:\n".format(tag))
        for rid in id_list:
            if terminalize:
                ntrw = taxalotl_config.get_resource_by_id(rid)
                if ntrw.id not in written:
                    ntrw.write_status(out_stream, indent="")
                written.add(ntrw.id)
                trw = taxalotl_config.get_terminalized_res_by_id(rid, "")
                print(type(trw))
                if trw is not ntrw and trw.id not in written:
                    trw.write_status(out_stream, indent="  ")
                    written.add(trw.id)
            else:
                rw = taxalotl_config.get_resource_by_id(rid)
                indent = "  " if rw.base_id in par_id_set else ""
                # out_stream.write('\n\nrid={}\n'.format(rid))
                if rw.id not in written:
                    rw.write_status(out_stream, indent=indent)
                    written.add(rw.id)
                if rw.is_abstract:
                    par_id_set.add(rw.id)


def unpack_resources(taxalotl_config, id_list):
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, "unpack")
        if not rw.has_been_downloaded():
            m = "{} will be downloaded first..."
            _LOG.info(m.format(rw.id))
            download_resources(taxalotl_config, [rw.id])
        if rw.has_been_unpacked():
            m = "{} was already present at {}"
            _LOG.info(m.format(rw.id, rw.unpacked_filepath))
        else:
            rw.unpack()


def normalize_resources(taxalotl_config, id_list):
    for rid in id_list:
        with VirtCommand(name="analyze-update", res_id=rid):
            rw = taxalotl_config.get_terminalized_res_by_id(rid, "normalize")
            if not rw.has_been_unpacked():
                m = "{} will be unpacked first..."
                _LOG.info(m.format(rw.id))
                unpack_resources(taxalotl_config, [rw.id])
            if rw.has_been_normalized():
                m = "{} was already normalized at {}"
                _LOG.info(m.format(rw.id, rw.normalized_filedir))
            else:
                rw.normalize()


def _iter_norm_term_res_internal_level_pairs(
    taxalotl_config, id_list, level_list, cmd_name
):
    """iterates over (non abstract resource, level) pairs

    Several cmds work on normalized resources and work on levels.
    This generator serves as a common iterator for them.
    Working on the specified and (as the inner loop) over the requested levels.
    """
    if level_list == [None]:
        level_list = PREORDER_PART_LIST
    for rid in id_list:
        res = taxalotl_config.get_terminalized_res_by_id(rid, cmd_name)
        if not res.has_been_normalized():
            normalize_resources(taxalotl_config, [rid])
        for part_name_to_split in level_list:
            if not NAME_TO_PARTS_SUBSETS[part_name_to_split]:
                _LOG.info(
                    '"{}" is a terminal group in the primary partition map'.format(
                        part_name_to_split
                    )
                )
            else:
                yield res, part_name_to_split


def info_on_resources(taxalotl_config, id_list, level_list):
    for res, part_name_to_split in _iter_norm_term_res_internal_level_pairs(
        taxalotl_config, id_list, level_list, "partition"
    ):
        write_info_for_res(out_stream, res, part_name_to_split)


def partition_resources(taxalotl_config, id_list, level_list):
    for res, part_name_to_split in _iter_norm_term_res_internal_level_pairs(
        taxalotl_config, id_list, level_list, "partition"
    ):
        with VirtCommand("partition", res_id=res.id, level=part_name_to_split):
            with use_tax_partitions():
                do_partition(res, part_name_to_split)


def exec_or_runtime_error(invocation, working_dir="."):
    rc = subprocess.call(invocation, cwd=working_dir)
    if rc != 0:
        qi = '", "'.join(invocation)
        m = 'Command\n"{}"\nfailed with returncode={}\n'
        raise RuntimeError(m.format(qi, rc))


def clone_otifacts(otifacts_dir):
    otifacts_url = "git@github.com:mtholder/OTifacts.git"
    m = 'Expecting OTifacts to be cloned at "{}". Will attempt to clone it from {}...'
    _LOG.warning(m.format(otifacts_dir, otifacts_url))
    exec_or_runtime_error(["git", "clone", otifacts_url, otifacts_dir])


def git_pull_otifacts(otifacts_dir):
    exec_or_runtime_error(["git", "pull"], working_dir=otifacts_dir)


def pull_otifacts(taxalotl_config):
    dest_dir = taxalotl_config.resources_dir
    taxalotl_dir = os.path.split(os.path.abspath(dest_dir))[0]
    repo_dir = os.path.split(taxalotl_dir)[0]
    otifacts_dir = os.path.join(repo_dir, "OTifacts")
    if not os.path.isdir(otifacts_dir):
        _LOG.debug(f"cloning to {otifacts_dir}")
        clone_otifacts(otifacts_dir)
    else:
        _LOG.debug(f"pulling to refresh {otifacts_dir}")
        git_pull_otifacts(otifacts_dir)
    all_res = read_all_otifacts(otifacts_dir)
    for res_type in [
        "external taxonomy",
        "open tree taxonomy",
        "id list",
        "open tree taxonomy idlist",
        "id to ncbi mapping",
    ]:
        ext_tax = filter_otifacts_by_type(all_res, res_type)
        by_root_id = partition_otifacts_by_root_element(ext_tax)
        for root_key, res_dict in by_root_id.items():
            fp = os.path.join(dest_dir, root_key + ".json")
            with OutFile(fp) as outs:
                write_as_json(res_dict, outs, indent=2)


def accumulate_taxon_dir_names(top_dir, name_to_paths):
    for root, dirs, files in os.walk(top_dir):
        if root.endswith(MISC_DIRNAME):
            continue
        if INP_TAXONOMY_DIRNAME in dirs or MISC_DIRNAME in dirs:
            name = os.path.split(root)[-1]
            name_to_paths.setdefault(name, []).append(root)


def cache_separator_names(taxalotl_config):
    rw = taxalotl_config.get_terminalized_res_by_id("ott", "")
    n2p = {}
    accumulate_taxon_dir_names(rw.partitioned_filepath, n2p)
    xl = list(n2p.keys())
    xl.sort()
    outfn = os.path.join(rw.partitioned_filepath, SEP_NAMES)
    with OutFile(outfn) as outs:
        write_as_json(xl, outs)
    _LOG.info("Separator dir names written to {}".format(outfn))
    outfn = os.path.join(rw.partitioned_filepath, SEP_MAPPING)
    for k, v in n2p.items():
        if len(v) > 1:
            _LOG.info("separator {} has multiple dirs: {}".format(k, v))
    with OutFile(outfn) as outs:
        write_as_json(n2p, outs)
    _LOG.info("Separator name to dir mapping written to {}".format(outfn))


def clean_resources(taxalotl_config, action, id_list, levels=None):
    if levels is None:
        levels = [None]
    if not id_list:
        rw = taxalotl_config.get_terminalized_res_by_id("ott", None)
        fp = os.path.join(rw.partitioned_filepath, GEN_MAPPING_FILENAME)
        if action == "build-partition-maps":
            if os.path.exists(fp):
                unlink(fp)
            else:
                _LOG.info(
                    'Mapping file "{}" does not exist. Skipping clean step'.format(fp)
                )
        else:
            raise NotImplementedError("clean of {} not yet implemented".format(action))
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, "clean")
        if action == "partition":
            if rw.has_been_partitioned():
                _LOG.info("Cleaning partition artifact for {}...".format(rid))
                rw.remove_partition_artifacts()
            else:
                _LOG.info(
                    "{} had not been partitioned. Skipping clean step...".format(rid)
                )
        elif action == "normalize":
            if rw.has_been_normalized():
                _LOG.info("Cleaning normalize artifact for {}...".format(rid))
                rw.remove_normalize_artifacts()
            else:
                _LOG.info(
                    "{} had not been normalized. Skipping clean step...".format(rid)
                )
        else:
            raise NotImplementedError("clean of {} not yet implemented".format(action))
