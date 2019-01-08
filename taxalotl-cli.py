#!/usr/bin/env python
from __future__ import print_function

import io
import sys
import os

from peyotl import (get_logger, read_as_json)

from taxalotl import TaxalotlConfig
from taxalotl.commands import (accumulate_separated_descendants,
                               analyze_update,
                               build_partition_maps,
                               cache_separator_names,
                               clean_resources,
                               compare_taxonomies,
                               diagnose_new_separators,
                               enforce_new_separators,
                               download_resources,
                               info_on_resources,
                               normalize_resources,
                               partition_resources,
                               pull_otifacts,
                               status_of_resources,
                               unpack_resources,
                               SEP_NAMES,
                               )
from taxalotl.partitions import (PART_NAMES,
                                 NAME_TO_PARTS_SUBSETS,
                                 NONTERMINAL_PART_NAMES,
                                 TERMINAL_PART_NAMES, )

_LOG = get_logger(__name__)

# Commands that don't take a resource ID
res_indep_cmds = ['build-partition-maps',
                  'cache-separator-names',
                  'clean-separation',
                  'compare-taxonomies',
                  'diagnose-new-separators',
                  'pull-otifacts',
                  ]
# Commands that take any resource ID
res_dep_cmds = ['accumulate-separated-descendants',
                'check-partition',
                'download',
                'enforce-new-separators',
                'info',
                'normalize',
                'partition',
                'status',
                'unpack',
                ]
# Commands that take an resource ID for a class of input resource (no version number suffix).
ver_inp_res_dep_cmds = ['analyze-update', ]
all_cmds = res_dep_cmds + res_indep_cmds + ver_inp_res_dep_cmds


def _validate_level_arg(taxalotl_config, level):
    if (level is not None) and (level not in NAME_TO_PARTS_SUBSETS):
        sep_dict = taxalotl_config.get_separator_dict()
        if level not in sep_dict:
            sn = set(NAME_TO_PARTS_SUBSETS.keys()).union(sep_dict.keys())
            raise RuntimeError('--level should be one of "{}"'.format('", "'.join(sn)))
    return True

def main_post_parse(args):
    taxalotl_config = TaxalotlConfig(filepath=args.config)
    hist_file = os.path.expanduser("~/.taxalotl_history")
    if os.path.isfile(hist_file):
        with io.open(hist_file, 'a', encoding='utf-8') as hout:
            hout.write('"{}"\n'.format('" "'.join(sys.argv)))
    try:
        if args.which == 'analyze-update':
            analyze_update(taxalotl_config, args.resources)
        elif args.which == 'clean-partition':
            clean_resources(taxalotl_config, 'partition', args.resources)
        elif args.which == 'clean-separation':
            clean_resources(taxalotl_config, 'separation', [], [args.level])
        elif args.which == 'cache-separator-names':
            cache_separator_names(taxalotl_config)
        elif args.which == 'compare-taxonomies':
            compare_taxonomies(taxalotl_config, [args.level])
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
        elif args.which == 'accumulate-separated-descendants':
            accumulate_separated_descendants(taxalotl_config, args.resources)
        elif args.which == 'pull-otifacts':
            pull_otifacts(taxalotl_config)
        elif args.which == 'diagnose-new-separators':
            _validate_level_arg(taxalotl_config, args.level)
            diagnose_new_separators(taxalotl_config, [args.level])
        elif args.which == 'enforce-new-separators':
            _validate_level_arg(taxalotl_config, args.level)
            enforce_new_separators(taxalotl_config, args.resources, [args.level])
        elif args.which == 'build-partition-maps':
            build_partition_maps(taxalotl_config)
        elif args.which == 'partition':
            if args.level is not None and args.level not in NAME_TO_PARTS_SUBSETS:
                raise RuntimeError('--level should be one of "{}"'.format('", "'.join(PART_NAMES)))
            partition_resources(taxalotl_config, args.resources, [args.level])
        elif args.which == 'info':
            if args.level is not None and args.level not in NAME_TO_PARTS_SUBSETS:
                raise RuntimeError('--level should be one of "{}"'.format('", "'.join(PART_NAMES)))
            info_on_resources(taxalotl_config, args.resources, [args.level])
        elif args.which == 'all':
            m = 'Currently you must enter a command to run. Use the --help option or see the Tutorial.md\n'
            sys.stdout.write(m)
            return 1
        else:
            raise NotImplementedError('"{}" action not implemented yet'.format(args.which))
    except Exception as x:
        if taxalotl_config.crash_with_stacktraces:
            raise
        sys.exit('taxalotl-cli: Exiting with exception:\n{}'.format(x))
    return 0


def main():
    import argparse

    description = "The main CLI for taxalotl"
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--config", type=str, help="the taxalotl.conf filepath (optional)")
    p.add_argument("--show-completions",
                   action="store_true",
                   default=False,
                   help="print the list of options for the next word in the command line")

    p.set_defaults(which="all")
    subp = p.add_subparsers(help="command help")
    # ANALYZE UPDATE
    analyze_update_p = subp.add_parser('analyze-update',
                                       help="calculates a diff between the last version of a "
                                            "taxonomy used and the latest version downloaded.")
    analyze_update_p.add_argument('resources', nargs="*", help="IDs of the resources to analyzed.")
    analyze_update_p.set_defaults(which="analyze-update")

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
                          help="Report only on the terminalized resource of each type.")
    status_p.set_defaults(which="status")
    # CACHE-separator-names
    compare_tax_p = subp.add_parser('compare-taxonomies',
                                    help="compare taxonomies for a separated dir")
    compare_tax_p.add_argument("--level",
                               default=None,
                               required=True,
                               help="The level of the taxonomy to compare")
    compare_tax_p.set_defaults(which="compare-taxonomies")

    # CACHE-separator-names
    cache_p = subp.add_parser('cache-separator-names',
                              help="Accumulate a list of separator names for tab-completion")
    cache_p.set_defaults(which="cache-separator-names")
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
    partition_p.add_argument('resources', nargs="+", help="IDs of the resources to partitition")
    partition_p.add_argument("--level",
                             default=None,
                             help="The level of the taxonomy to partition")

    partition_p.set_defaults(which="partition")

    # INFO
    info_p = subp.add_parser('info',
                                  help="Report statistics about a resource")
    info_p.add_argument('resources', nargs="+", help="IDs of the resources")
    info_p.add_argument("--level",
                             default=None,
                             help="The level of the taxonomy to partition")

    info_p.set_defaults(which="info")

    # DIAGNOSE-NEW-SEPARATORS
    diag_sep_p = subp.add_parser('diagnose-new-separators',
                                 help="Uses the last OTT build to find taxa IDs that "
                                      "feature are common to the relevant inputs")
    diag_sep_p.add_argument("--level",
                            default=None,
                            help="The partition that is the root of the separation (default all)")

    diag_sep_p.set_defaults(which="diagnose-new-separators")
    # ENFORCE-NEW-SEPARATORS
    enf_sep_p = subp.add_parser('enforce-new-separators',
                                help="Uses the __sep__.json files created by "
                                     "diagnose-new-separators to partition by unproblematic "
                                     "taxa")
    enf_sep_p.add_argument('resources', nargs="*", help="IDs of the resources to separate")
    enf_sep_p.add_argument("--level",
                           default=None,
                           help="The partition that is the root of the separation (default all)")
    enf_sep_p.set_defaults(which="enforce-new-separators")
    # ACCUMULATE-SEPARATED-DESCENDANTS
    accum_sep_des_p = subp.add_parser('accumulate-separated-descendants',
                                      help="Should be run after enforce-separators and before compare-taxonomies")
    accum_sep_des_p.add_argument('resources', nargs="*", help="IDs of the resources")
    accum_sep_des_p.set_defaults(which="accumulate-separated-descendants")

    # BUILD-PARTITION-MAPS
    build_partition_maps_p = subp.add_parser('build-partition-maps',
                                             help="Uses the last OTT build to find the "
                                                  "ID mappings needed to "
                                                  "partition the inputs taxonomies.")
    build_partition_maps_p.set_defaults(which="build-partition-maps")
    # CLEAN-PARTITION
    clean_p = subp.add_parser('clean-partition',
                              help="remove the results the partition+enforce-new-separator for a resource.")
    clean_p.add_argument('resources', nargs="*", help="IDs of the resources to clean")
    clean_p.set_defaults(which='clean-partition')
    # CLEAN-PARTITION
    clean_s_p = subp.add_parser('clean-separation',
                              help="remove the results the diagnose-new-separator for a resource.")
    clean_s_p.add_argument("--level",
                           default=None,
                           help="The partition that is the root of the separation (default all)")

    clean_s_p.set_defaults(which='clean-separation')

    # Handle --show-completions differently from the others, because
    #   argparse does not help us out here... at all
    if "--show-completions" in sys.argv:
        a = sys.argv[1:]
        univ = frozenset(['--config', ])
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
            if sel_cmd in res_dep_cmds \
              or sel_cmd in ['compare-taxonomies'] \
              or sel_cmd in ver_inp_res_dep_cmds:
                # From Ned Batchelder's answer on http://stackoverflow.com/a/14728477
                class ArgumentParserError(Exception):
                    pass

                # noinspection PyClassHasNoInit
                class ThrowingArgumentParser(argparse.ArgumentParser):
                    def error(self, message):
                        raise ArgumentParserError(message)

                fake_parser = ThrowingArgumentParser()
                fake_parser.add_argument("--config", type=str)
                fake_parser.add_argument('blah', nargs="*")
                comp_list = []
                taxalotl_config = None
                try:
                    fa = fake_parser.parse_known_args()[0]
                    config = fa.config
                    taxalotl_config = TaxalotlConfig(filepath=config)
                    if sel_cmd in res_dep_cmds:
                        comp_list = list(taxalotl_config.resources_mgr.resources.keys())
                    elif sel_cmd in ver_inp_res_dep_cmds:
                        comp_list = list(
                            taxalotl_config.resources_mgr.abstract_input_resource_types())
                except Exception as _excep:
                    _LOG.warn('Exception: {}'.format(_excep))
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
                elif sel_cmd in ('diagnose-new-separators', 'enforce-new-separators'):
                    # sys.stderr.write(str(a))
                    if '--level' == a[-1] or (len(a) > 1 and '--level' == a[-2]):
                        comp_list = list(TERMINAL_PART_NAMES)
                    elif '--level' not in a:
                        comp_list.extend(['--level'])
                elif sel_cmd in ['compare-taxonomies']:
                    rw = taxalotl_config.get_terminalized_res_by_id("ott", '')
                    outfn = os.path.join(rw.partitioned_filepath, SEP_NAMES)
                    if os.path.exists(outfn):
                        comp_list.extend(read_as_json(outfn))

        sys.stdout.write('{}\n'.format(' '.join(comp_list)))
    else:
        rc = main_post_parse(p.parse_args())
        sys.exit(rc)


if __name__ == "__main__":
    main()
