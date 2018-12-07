#!/usr/bin/env python
from __future__ import print_function

import subprocess
import os
import sys

from peyotl import (get_logger,
                    read_all_otifacts,
                    read_as_json,
                    filter_otifacts_by_type,
                    partition_otifacts_by_root_element,
                    write_as_json)

from taxalotl.partitions import (GEN_MAPPING_FILENAME,
                                 check_partition,
                                 do_partition,
                                 NAME_TO_PARTS_SUBSETS,
                                 PART_NAMES,
                                 PREORDER_PART_LIST,
                                 TERMINAL_PART_NAMES)
from taxalotl.tax_partition import use_tax_partitions
from taxalotl.dynamic_partitioning import perform_dynamic_separation
from taxalotl.compare import compare_taxonomies_in_dir

_LOG = get_logger(__name__)
out_stream = sys.stdout

SEP_NAMES = '__separator_names__.json'
SEP_MAPPING = '__separator_names_to_dir__.json'

def analyze_update_to_resources(prev, curr):
    m = 'Could not analyze taxonomy update because {} has not been partitioned.'
    for el in [prev, curr]:
        if not el.has_been_partitioned():
            raise RuntimeError(m.format(el.id))
    raise NotImplementedError('ha')

def analyze_update(taxalotl_config, id_list):
    ott_res = taxalotl_config.get_terminalized_res_by_id('ott')
    ott_inps = ott_res.input_dict()
    for i in id_list:
        prev_used_version_id = ott_inps.get(i)
        if prev_used_version_id is None:
            _LOG.info('Taxonomy {} was a not an input to the current version of OTT, so an update cannot be analyzed.'.format(i))
            continue
        x = taxalotl_config.get_all_terminalized_res_by_id(i)
        xi = [j.id for j in x]
        try:
            prev_used_version_index = xi.index(prev_used_version_id)
        except ValueError:
            _LOG.info('The previously used version of {} is {}, but this is not known to taxalotl, you might need to run the pull-otifacts action.'.format(i, prev_used_version_id))
            continue
        if prev_used_version_index == len(xi) - 1:
            _LOG.info('The previously used version of {} is {} appears to be the latest version of that resource, you might need to run the pull-otifacts action.'.format(i, prev_used_version_id))
            continue
        prev_used = x[prev_used_version_index]
        latest_wrapper = x[-1]
        analyze_update_to_resources(prev_used, latest_wrapper)

def download_resources(taxalotl_config, id_list):
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, 'download')
        lic_urls, lic_tou = taxalotl_config.get_known_license_info(rw)
        if lic_urls or lic_tou:
            prompt = 'The following license-related URLs were found:\n  '
            prompt += '\n  '.join(lic_urls)
            prompt += '\n\nThe following license or terms of use information was found:\n\n  '
            prompt += '\n\n  '.join(lic_tou)
        else:
            prompt = 'No stored license info based on the last pull from the OTifacts repository.'
        tag = '\n\nOther license or terms of use may apply if the OTifacts repository is not ' \
              'up-to-date.\nEnter y to continue downloading {} : '.format(rid)
        prompt += tag
        _LOG.info(repr(prompt))
        try:
            # noinspection PyCompatibility
            resp = input(prompt)
        except NameError:
            resp = input(prompt)
        if resp != 'y':
            _LOG.info('download of {} skipped due to lack of affirmative response.'.format(rid))
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
    return [["abstract classes", a_list],
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

def status_of_resources(taxalotl_config,
                        id_list,
                        ids_only=False,
                        by_status=False,
                        terminal_only=False):
    terminalize = True
    if not id_list:
        id_list = get_list_of_all_resources(taxalotl_config)
        if not id_list:
            return
        terminalize = False
    if terminal_only:
        x = []
        for i in id_list:
            ri = taxalotl_config.get_terminalized_res_by_id(i, '')
            if ri.id == i:
                x.append(i)
        id_list = x
    res = taxalotl_config.resources_mgr.resources
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
    par_id_set = set()
    for tag, id_list in t_and_id_list:
        if tag:
            out_stream.write("{}:\n".format(tag))
        for rid in id_list:
            if terminalize:
                ntrw = taxalotl_config.get_resource_by_id(rid)
                ntrw.write_status(out_stream, indent='')
                trw = taxalotl_config.get_terminalized_res_by_id(rid, '')
                trw.write_status(out_stream, indent='  ')
            else:
                rw = taxalotl_config.get_resource_by_id(rid)
                indent = '  ' if rw.base_id in par_id_set else ''
                rw.write_status(out_stream, indent=indent)
                if rw.is_abstract:
                    par_id_set.add(rw.id)


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
        res = taxalotl_config.get_terminalized_res_by_id(rid, 'partition')
        if not res.has_been_normalized():
            normalize_resources(taxalotl_config, [rid])
        for part_name_to_split in level_list:
            if not NAME_TO_PARTS_SUBSETS[part_name_to_split]:
                _LOG.info('"{}" is a terminal group in the primary partition map'.format(
                    part_name_to_split))
                continue
            with use_tax_partitions() as tax_cache:
                do_partition(res, part_name_to_split)


def check_partition_resources(taxalotl_config, id_list, level_list):
    if level_list == [None]:
        level_list = PREORDER_PART_LIST
    for rid in id_list:
        res = taxalotl_config.get_terminalized_res_by_id(rid, 'partition')
        if not res.has_been_normalized():
            normalize_resources(taxalotl_config, [rid])
        for part_name_to_split in level_list:
            if not NAME_TO_PARTS_SUBSETS[part_name_to_split]:
                _LOG.info('"{}" is a terminal group in the primary partition map'.format(
                    part_name_to_split))
                continue
            with use_tax_partitions() as tax_cache:
                check_partition(res, part_name_to_split)

def exec_or_runtime_error(invocation, working_dir='.'):
    rc = subprocess.call(invocation, cwd=working_dir)
    if rc != 0:
        qi = '", "'.join(invocation)
        m = 'Command\n"{}"\nfailed with returncode={}\n'
        raise RuntimeError(m.format(qi, rc))

def clone_otifacts(otifacts_dir):
    otifacts_url = 'git@github.com:mtholder/OTifacts.git'
    m = 'Expecting OTifacts to be cloned at "{}". Will attempt to clone it from {}...'
    _LOG.warn(m.format(otifacts_dir, otifacts_url))
    exec_or_runtime_error(['git', 'clone', otifacts_url, otifacts_dir])


def git_pull_otifacts(otifacts_dir):
    exec_or_runtime_error(['git', 'pull'], working_dir=otifacts_dir)

def pull_otifacts(taxalotl_config):
    dest_dir = taxalotl_config.resources_dir
    taxalotl_dir = os.path.split(os.path.abspath(dest_dir))[0]
    repo_dir = os.path.split(taxalotl_dir)[0]
    otifacts_dir = os.path.join(repo_dir, 'OTifacts')
    if not os.path.isdir(otifacts_dir):
        clone_otifacts(otifacts_dir)
    else:
        git_pull_otifacts(otifacts_dir)
    all_res = read_all_otifacts(otifacts_dir)
    for res_type in ['external taxonomy', 'open tree taxonomy', 'id list']:
        ext_tax = filter_otifacts_by_type(all_res, res_type)
        by_root_id = partition_otifacts_by_root_element(ext_tax)
        for root_key, res_dict in by_root_id.items():
            fp = os.path.join(dest_dir, root_key + '.json')
            write_as_json(res_dict, fp, indent=2, separators=(',', ': '))


NEW_SEP_FILENAME = '__sep__.json'


def diagnose_new_separators(taxalotl_config, level_list):
    rw = taxalotl_config.get_terminalized_res_by_id("ott", 'diagnose-new-separators')
    if not rw.has_been_partitioned():
        partition_resources(taxalotl_config, ["ott"], PREORDER_PART_LIST)
    pd = rw.partitioned_filepath
    if level_list == [None]:
        level_list = PART_NAMES
    for part_name in level_list:
        nsd = rw.diagnose_new_separators(current_partition_key=part_name)
        if not nsd:
            _LOG.info("no new separtors in {}.".format(part_name))
        else:
            for k, sd in nsd.items():
                _LOG.info('{} new separators in {}'.format(sd.num_separators(),
                                                           part_name))
                fp = os.path.join(pd, k, NEW_SEP_FILENAME)
                write_as_json(sd.as_dict(), fp, sort_keys=True, indent=2)
                _LOG.info("new separators written to {}".format(fp))


def enforce_new_separators(taxalotl_config, id_list, level_list):
    if level_list == [None]:
        level_list = list(PREORDER_PART_LIST) + list(TERMINAL_PART_NAMES)
    ott_res = taxalotl_config.get_terminalized_res_by_id("ott", 'enforce-new-separators')
    if not ott_res.has_been_partitioned():
        partition_resources(taxalotl_config, ["ott"], PREORDER_PART_LIST)
    if not id_list:
        id_list = ["ott", "gbif", "irmng", "silva", "worms", "ncbi"]
    for src in id_list:
        if src is not None:
            res = taxalotl_config.get_terminalized_res_by_id(src, 'enforce-new-separators')
            if not res.has_been_partitioned():
                partition_resources(taxalotl_config, [src], PREORDER_PART_LIST)
        with use_tax_partitions() as cache:
            for part_name in level_list:
                perform_dynamic_separation(ott_res,
                                           src_id=src,
                                           part_key=part_name,
                                           sep_fn=NEW_SEP_FILENAME,
                                           suppress_cache_flush=True)


def build_partition_maps(taxalotl_config):
    rw = taxalotl_config.get_terminalized_res_by_id("ott", 'partition')
    if not rw.has_been_partitioned():
        partition_resources(taxalotl_config, ["ott"], PREORDER_PART_LIST)
    nsd = rw.build_paritition_maps()
    if not nsd:
        return
    pd = rw.partitioned_filepath
    mfp = os.path.join(pd, GEN_MAPPING_FILENAME)
    write_as_json(nsd, mfp, indent=2)
    _LOG.info("Partitions maps written to {}".format(mfp))


def accumulate_taxon_dir_names(top_dir, name_to_paths):
    for root, dirs, files in os.walk(top_dir):
        if root.endswith('__misc__'):
            continue
        if '__inputs__' in dirs or '__misc__' in dirs:
            name = os.path.split(root)[-1]
            name_to_paths.setdefault(name, []).append(root)


def cache_separator_names(taxalotl_config):
    rw = taxalotl_config.get_terminalized_res_by_id("ott", '')
    n2p = {}
    accumulate_taxon_dir_names(rw.partitioned_filepath, n2p)
    xl = list(n2p.keys())
    xl.sort()
    outfn = os.path.join(rw.partitioned_filepath, SEP_NAMES)
    write_as_json(xl, outfn)
    _LOG.info("Separator dir names written to {}".format(outfn))
    outfn = os.path.join(rw.partitioned_filepath, SEP_MAPPING)
    for k, v in n2p.items():
        if len(v) > 1:
            _LOG.info("separator {} has multiple dirs: {}".format(k, v))
    write_as_json(n2p, outfn)
    _LOG.info("Separator name to dir mapping written to {}".format(outfn))


def compare_taxonomies(taxalotl_config, levels):
    assert levels != [None]
    rw = taxalotl_config.get_terminalized_res_by_id("ott", '')
    outfn = os.path.join(rw.partitioned_filepath, SEP_MAPPING)
    if not os.path.exists(outfn):
        cache_separator_names(taxalotl_config)
    todir = read_as_json(outfn)
    for level in levels:
        try:
            tax_dir_list = todir[level]
        except KeyError:
            raise ValueError('The level "{}" is not separator name'.format(level))
        for tax_dir in tax_dir_list:
            m = 'Will compare taxonomies for "{}" based on {}'
            _LOG.info(m.format(level, tax_dir))
            compare_taxonomies_in_dir(taxalotl_config, tax_dir)


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


def accumulate_separated_descendants(taxalotl_config, id_list):
    rw = taxalotl_config.get_terminalized_res_by_id("ott", '')
    outfn = os.path.join(rw.partitioned_filepath, SEP_MAPPING)
    if not os.path.exists(outfn):
        cache_separator_names(taxalotl_config)
    assert os.path.exists(outfn)
    part_to_dir = read_as_json(outfn)
    # created a post-order list by using the length of the directory...
    dir_tuple_list = []
    for d in part_to_dir.values():
        for ds in d:
            dir_tuple_list.append((len(ds), ds))
    dir_tuple_list.sort(reverse=True)
    postorder = [i[1] for i in dir_tuple_list]
    for i in id_list:
        _LOG.info('accumulate_separated_descendants for {}'.format(i))
        res = taxalotl_config.get_terminalized_res_by_id(i, '')
        for d in postorder:
            _LOG.info('accumulate_separated_descendants for {}'.format(d))
            res.accumulate_separated_descendants(d)
