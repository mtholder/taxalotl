#!/usr/bin/env python
from __future__ import print_function

import os
import sys

from peyotl import (get_logger,
                    read_all_otifacts,
                    filter_otifacts_by_type,
                    partition_otifacts_by_root_element,
                    write_as_json)

from taxalotl import TaxalotlConfig
from taxalotl.partitions import (GEN_MAPPING_FILENAME,
                                 PART_NAMES,
                                 PART_FRAG_BY_NAME,
                                 PARTS_BY_NAME,
                                 NONTERMINAL_PART_NAMES,
                                 PREORDER_PART_LIST)

_LOG = get_logger(__name__)
out_stream = sys.stdout
res_indep_cmds = ['pull-otifacts', 'diagnose-new-separators', 'build-partition-maps']
res_dep_cmds = ['clean', 'status', 'download', 'unpack', 'normalize', 'partition']
all_cmds = res_dep_cmds + res_indep_cmds


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
    for res_type in ['external taxonomy', 'open tree taxonomy']:
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


def main_post_parse(args):
    taxalotl_config = TaxalotlConfig(filepath=args.config, resources_dir=args.resources_dir)
    try:
        if args.which == 'clean':
            if args.action not in all_cmds:
                m = "Expecting clean action to be one of: {}"
                raise ValueError(m.format(', '.join(all_cmds)))
            clean_resources(taxalotl_config, args.action, args.resources)
        elif args.which == 'download':
            download_resources(taxalotl_config, args.resources)
        elif args.which == 'status':
            status_of_resources(taxalotl_config,
                                args.resources,
                                ids_only=args.ids_only,
                                by_status=args.by_status,
                                terminal_only=args.terminal)
        elif args.which == 'unpack':
            unpack_resources(taxalotl_config, args.resources)
        elif args.which == 'normalize':
            normalize_resources(taxalotl_config, args.resources)
        elif args.which == 'pull-otifacts':
            pull_otifacts(taxalotl_config)
        elif args.which == 'diagnose-new-separators':
            diagnose_new_separators(taxalotl_config)
        elif args.which == 'build-partition-maps':
            build_partition_maps(taxalotl_config)
        elif args.which == 'partition':
            if args.level is not None and args.level not in PARTS_BY_NAME:
                raise RuntimeError('--level should be one of "{}"'.format('", "'.join(PART_NAMES)))
            partition_resources(taxalotl_config, args.resources, [args.level])
        else:
            raise NotImplementedError('"{}" action not implemented yet'.format(args.which))
    except Exception as x:
        if taxalotl_config.crash_with_stacktraces:
            raise
        sys.exit('taxalotl-cli: Exiting with exception:\n{}'.format(x))


def main():
    import argparse

    description = "The main CLI for taxalotl"
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--resources-dir", type=str, help="the resources directory (optional)")
    p.add_argument("--config", type=str, help="the taxalotl.conf filepath (optional)")
    p.add_argument("--show-completions",
                   action="store_true",
                   default=False,
                   help="print the list of options for the next word in the command line")

    p.set_defaults(which="all")
    subp = p.add_subparsers(help="command help")
    # PULL OTifacts
    pull_otifacts_p = subp.add_parser('pull-otifacts',
                                      help="refresh list of taxonomic artifacts from OTifacts repo")
    pull_otifacts_p.set_defaults(which="pull-otifacts")
    # STATUS
    status_p = subp.add_parser('status',
                               help="report the status of a resource (or all resources)")
    status_p.add_argument('resources', nargs="*", help="IDs of the resources to report status on")
    status_p.add_argument("-i", "--ids-only",
                          action='store_true',
                          default=False,
                          help="just list the IDs")
    status_p.add_argument("--by-status",
                          action='store_true',
                          default=False,
                          help="group the report by status")
    status_p.add_argument("--terminal",
                          action='store_true',
                          default=False,
                          help="Only report on the most recent, terminalized resource of each type.")
    status_p.set_defaults(which="status")
    # DOWNLOAD
    download_p = subp.add_parser('download', help="download an artifact to your local filesystem")
    download_p.add_argument('resources', nargs="+", help="IDs of the resources to download")
    download_p.set_defaults(which="download")
    # UNPACK
    unpack_p = subp.add_parser('unpack',
                               help="unpack an resource (downloads if necessary)")
    unpack_p.add_argument('resources', nargs="+", help="IDs of the resources to unpack")
    unpack_p.set_defaults(which="unpack")
    # NORMALIZE
    normalize_p = subp.add_parser('normalize',
                                  help="converts to the OTT format (unpacks if necessary)")
    normalize_p.add_argument('resources', nargs="+", help="IDs of the resources to normalize")
    normalize_p.set_defaults(which="normalize")
    # PARTITION
    partition_p = subp.add_parser('partition',
                                  help="Breaks the resource taxon")
    partition_p.add_argument('resources', nargs="+", help="IDs of the resources to unpack")
    partition_p.add_argument("--level",
                             default=None,
                             help="The level of the taxonomy to partition")

    partition_p.set_defaults(which="partition")
    # DIAGNOSE-NEW-SEPARATORS
    diag_sep_p = subp.add_parser('diagnose-new-separators',
                                 help="Uses the last OTT build to find taxa IDs that "
                                      "feature are common to the relevant inputs")
    diag_sep_p.set_defaults(which="diagnose-new-separators")
    # BUILD-PARTITION-MAPS
    build_partition_maps_p = subp.add_parser('build-partition-maps',
                                             help="Uses the last OTT build to find the "
                                                  "ID mappings needed to "
                                                  "partition the inputs taxonomies.")
    build_partition_maps_p.set_defaults(which="build-partition-maps")
    # CLEAN
    clean_p = subp.add_parser('clean', help="remove the results of an action for a resource.")
    clean_p.add_argument('action', help="command to clean up for")
    clean_p.add_argument('resources', nargs="*", help="IDs of the resources to clean")
    clean_p.set_defaults(which='clean')
    # Handle --show-completions differently from the others, because
    #   argparse does not help us out here... at all
    if "--show-completions" in sys.argv:
        a = sys.argv[1:]
        univ = frozenset(['--resources-dir', '--config'])
        sel_cmd = None
        num_cmds = 0
        for c in all_cmds:
            if c in a:
                if sel_cmd is None:
                    sel_cmd = c
                num_cmds += 1
        comp_list = []
        if sel_cmd is None:
            comp_list = []
            for u in univ:
                found = False
                for arg in a:
                    if arg.startswith(u):
                        found = True
                        break
                if not found:
                    comp_list.append(u)
            comp_list.extend(all_cmds)
        else:
            if sel_cmd in res_dep_cmds:
                # From Ned Batchelder's answer on http://stackoverflow.com/a/14728477
                class ArgumentParserError(Exception):
                    pass

                # noinspection PyClassHasNoInit
                class ThrowingArgumentParser(argparse.ArgumentParser):
                    def error(self, message):
                        raise ArgumentParserError(message)

                fake_parser = ThrowingArgumentParser()
                fake_parser.add_argument("--resources-dir", type=str)
                fake_parser.add_argument("--config", type=str)
                fake_parser.add_argument('blah', nargs="*")
                comp_list = []
                try:
                    fa, unk = fake_parser.parse_known_args()
                    resdir, config = fa.resources_dir, fa.config
                    taxalotl_config = TaxalotlConfig(filepath=config, resources_dir=resdir)
                    comp_list = list(taxalotl_config.resources_mgr.resources.keys())
                except:
                    pass

                if sel_cmd == 'status':
                    if '-i' not in a and '--ids-only' not in a:
                        comp_list.extend(["-i", "--ids-only"])
                    for x in ['--by-status', '--terminal']:
                        if x not in a:
                            comp_list.extend([x])
                elif sel_cmd == 'partition':
                    # sys.stderr.write(str(a))
                    if '--level' == a[-1] or (len(a) > 1 and '--level' == a[-2]):
                        comp_list = list(NONTERMINAL_PART_NAMES)
                    elif '--level' not in a:
                        comp_list.extend(['--level'])
                elif sel_cmd == 'clean' and num_cmds == 1:
                    all_cmds.remove('clean')
                    all_cmds.remove('status')
                    comp_list.extend(all_cmds)

        sys.stdout.write('{}\n'.format(' '.join(comp_list)))
    else:
        main_post_parse(p.parse_args())


if __name__ == "__main__":
    main()
