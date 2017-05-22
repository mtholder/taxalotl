#!/usr/bin/env python
from __future__ import print_function

import os
import sys

from peyotl import (get_logger,
                    read_all_otifacts,
                    filter_otifacts_by_type,
                    partition_otifacts_by_root_element,
                    write_as_json)

from taxalotl.partitions import (GEN_MAPPING_FILENAME,
                                 PART_FRAG_BY_NAME,
                                 PARTS_BY_NAME,
                                 PREORDER_PART_LIST)

_LOG = get_logger(__name__)
out_stream = sys.stdout


def download_resources(taxalotl_config, id_list):
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, 'download')
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
    return [["abstract classes", a_list],
            ["not downloaded", nd_list],
            ["downloaded, but not unpacked", dnu_list],
            ["unpacked, but not normalized", unn_list],
            ["normalized", n_list],
            ["parititioned", p_list],
            ]


def status_of_resources(taxalotl_config,
                        id_list,
                        ids_only=False,
                        by_status=False,
                        terminal_only=False):
    terminalize = True
    res = taxalotl_config.resources_mgr.resources
    if not id_list:
        terminalize = False
        id_list = list(res.keys())
        id_list.sort()
    if terminal_only:
        x = []
        for i in id_list:
            ri = taxalotl_config.get_terminalized_res_by_id(i, '')
            if ri.id == i:
                x.append(i)
        id_list = x
    if by_status:
        t_and_id_list = _group_by_status(res, id_list)
    else:
        t_and_id_list = [["", id_list]]
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
    for tag, id_list in t_and_id_list:
        if tag:
            out_stream.write("{}:\n".format(tag))
        for rid in id_list:
            if terminalize:
                rw = taxalotl_config.get_terminalized_res_by_id(rid, '')
            else:
                rw = taxalotl_config.get_resource_by_id(rid)
            rw.write_status(out_stream, indent=' ' * 4)


def unpack_resources(taxalotl_config, id_list):
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, 'unpack')
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
        rw = taxalotl_config.get_terminalized_res_by_id(rid, 'normalize')
        if not rw.has_been_unpacked():
            m = "{} will be unpacked first..."
            _LOG.info(m.format(rw.id))
            unpack_resources(taxalotl_config, [rw.id])
        if rw.has_been_normalized():
            m = "{} was already normalized at {}"
            _LOG.info(m.format(rw.id, rw.normalized_filepath))
        else:
            rw.normalize()


def partition_resources(taxalotl_config, id_list, level_list):
    if level_list == [None]:
        level_list = PREORDER_PART_LIST
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, 'partition')
        if not rw.has_been_unpacked():
            m = "{} will be unpacked first..."
            _LOG.info(m.format(rw.id))
            unpack_resources(taxalotl_config, [rw.id])
        for level in level_list:
            rw.partition(level, PARTS_BY_NAME[level], PART_FRAG_BY_NAME[level])


def pull_otifacts(taxalotl_config):
    dest_dir = taxalotl_config.resources_dir
    taxalotl_dir = os.path.split(os.path.abspath(dest_dir))[0]
    repo_dir = os.path.split(taxalotl_dir)[0]
    otifacts_dir = os.path.join(repo_dir, 'OTifacts')
    if not os.path.isdir(otifacts_dir):
        m = 'Expecting OTifacts to be cloned as sibling of this directory at "{}"'
        raise RuntimeError(m.format(otifacts_dir))
    all_res = read_all_otifacts(otifacts_dir)
    for res_type in ['external taxonomy', 'open tree taxonomy', 'id list']:
        ext_tax = filter_otifacts_by_type(all_res, res_type)
        by_root_id = partition_otifacts_by_root_element(ext_tax)
        for root_key, res_dict in by_root_id.items():
            fp = os.path.join(dest_dir, root_key + '.json')
            write_as_json(res_dict, fp, indent=2, separators=(',', ': '))


def diagnose_new_separators(taxalotl_config):
    partition_resources(taxalotl_config, ["ott"], PREORDER_PART_LIST)
    rw = taxalotl_config.get_terminalized_res_by_id("ott", 'partition')
    nsd = rw.diagnose_new_separators(current_partition_key='Glaucophyta')
    if not nsd:
        return
    for k, sd in nsd.items():
        _LOG.info('New separators {} => {}'.format(k, sd))


def build_partition_maps(taxalotl_config):
    partition_resources(taxalotl_config, ["ott"], PREORDER_PART_LIST)
    rw = taxalotl_config.get_terminalized_res_by_id("ott", 'partition')
    nsd = rw.build_paritition_maps()
    if not nsd:
        return
    pd = rw.partitioned_filepath
    mfp = os.path.join(pd, GEN_MAPPING_FILENAME)
    write_as_json(nsd, mfp, indent=2)
    _LOG.info("Partitions maps written to {}".format(mfp))


def clean_resources(taxalotl_config, action, id_list):
    if not id_list:
        if action == 'build-partition-maps':
            rw = taxalotl_config.get_terminalized_res_by_id("ott", None)
            fp = os.path.join(rw.partitioned_filepath, GEN_MAPPING_FILENAME)
            if os.path.exists(fp):
                _LOG.info('Removing "{}"...'.format(fp))
                os.unlink(fp)
            else:
                _LOG.info('Mapping file "{}" does not exist. Skipping clean step'.format(fp))
        else:
            raise NotImplementedError("clean of {} not yet implemented".format(action))
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, 'clean')
        if action == 'partition':
            if rw.has_been_partitioned():
                _LOG.info("Cleaning partition artifact for {}...".format(rid))
                rw.remove_partition_artifacts()
            else:
                _LOG.info("{} had not been partitioned. Skipping clean step...".format(rid))
        elif action == 'normalize':
            if rw.has_been_normalized():
                _LOG.info("Cleaning normalize artifact for {}...".format(rid))
                rw.remove_normalize_artifacts()
            else:
                _LOG.info("{} had not been normalized. Skipping clean step...".format(rid))
        else:
            raise NotImplementedError("clean of {} not yet implemented".format(action))
