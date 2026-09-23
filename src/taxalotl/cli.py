#!/usr/bin/env python
from __future__ import print_function

import os
import sys

from peyutil import read_as_json
import argparse
from . import TaxalotlConfig
from .commands import (
    # analyze_update,
    add_mapping,
    clean_resources,
    download_resources,
    grep_in_res,
    info_on_resources,
    normalize_resources,
    partition_resources,
    pull_otifacts,
    status_of_resources,
    unpack_resources,
)
from .cmds.partitions import (
    PART_NAMES,
    NAME_TO_PARTS_SUBSETS,
    NONTERMINAL_PART_NAMES,
    TERMINAL_PART_NAMES,
)
import logging

LOGLEVEL = os.environ.get("LOGLEVEL", "WARNING").upper()
logging.basicConfig(level=LOGLEVEL)
_LOG = logging.getLogger(__name__)

# Commands that don't take a resource ID
res_indep_cmds = [
    "compare-taxonomies",
    "pull-otifacts",
]
# Commands that take any resource ID
res_dep_cmds = [
    "add-mapping",
    "check-partition",
    "clean-partition",
    "download",
    "grep",
    "info",
    "normalize",
    "partition",
    "status",
    "unpack",
]
# Commands that take an resource ID for a class of input resource (no version number suffix).
ver_inp_res_dep_cmds = []
all_cmds = res_dep_cmds + res_indep_cmds + ver_inp_res_dep_cmds


def _verify_level_arg(lev_arg):
    if lev_arg is not None and lev_arg not in NAME_TO_PARTS_SUBSETS:
        opts = '", "'.join(PART_NAMES)
        raise RuntimeError(f'--level should be one of "{opts}"')
    return [lev_arg]


def main_post_parse(args):
    cfg = TaxalotlConfig(filepath=args.config)
    try:
        # if args.which == 'analyze-update':
        #     analyze_update(cfg, args.resources, [args.level])
        # elif
        if args.which == "clean-partition":
            clean_resources(cfg, "partition", args.resources)
        elif args.which == "download":
            download_resources(cfg, args.resources)
        elif args.which == "status":
            status_of_resources(
                cfg,
                args.resources,
                ids_only=args.ids_only,
                by_status=args.by_status,
                terminal_only=args.terminal,
            )
        elif args.which == "unpack":
            unpack_resources(cfg, args.resources)
        elif args.which == "normalize":
            normalize_resources(cfg, args.resources)
        elif args.which == "pull-otifacts":
            pull_otifacts(cfg)
        elif args.which == "partition":
            lev = _verify_level_arg(args.level)
            partition_resources(cfg, args.strategy, args.resources, lev)
        elif args.which == "info":
            lev = _verify_level_arg(args.level)
            info_on_resources(cfg, args.resources, lev)
        elif args.which == "grep":
            if args.name:
                if len(args.name) > 1:
                    raise RuntimeError("Only 1 name argument allowed")
                name_arg = args.name[0]
                if args.tax_id:
                    raise RuntimeError("name or tax_id_field can be used, not both")
                tax_id_arg = None
            elif args.tax_id:
                if len(args.tax_id) > 1:
                    raise RuntimeError("Only 1 tax_id argument allowed")
                tax_id_arg = args.tax_id[0]
                name_arg = None
            else:
                raise RuntimeError("either name or tax_id_field must be used.")
            grep_in_res(cfg, args.resources, name_arg, tax_id_arg, args.target)
        elif args.which == "add-mapping":
            add_mapping(cfg, args.ott_id, args.external_id)
        elif args.which == "all":
            m = "Currently you must enter a command to run. Use the --help option or see the Tutorial.md\n"
            sys.stdout.write(m)
            return 1
        else:
            raise NotImplementedError(
                '"{}" action not implemented yet'.format(args.which)
            )
    except Exception as x:
        if cfg.crash_with_stacktraces:
            raise
        sys.exit("taxalotl-cli: Exiting with exception:\n{}".format(x))
    return 0


def _add_level_arg(parser, req=False):
    parser.add_argument(
        "--level", default=None, required=req, help="The highest taxon to work on."
    )


def main():

    description = "The main CLI for taxalotl"
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--config", type=str, help="the taxalotl.conf filepath (optional)")
    p.add_argument(
        "--show-completions",
        action="store_true",
        default=False,
        help="print the list of options for the next word in the command line",
    )

    p.set_defaults(which="all")
    subp = p.add_subparsers(help="command help")
    # ANALYZE UPDATE
    # analyze_update_p = subp.add_parser('analyze-update',
    #                                    help="calculates a diff between the last version of a "
    #                                         "taxonomy used and the latest version downloaded.")
    # analyze_update_p.add_argument('resources', nargs=2, help="IDs of the resources to analyzed.")
    # _add_level_arg(analyze_update_p)
    # analyze_update_p.set_defaults(which="analyze-update")

    # PULL OTifacts
    pull_otifacts_p = subp.add_parser(
        "pull-otifacts", help="refresh list of taxonomic artifacts from OTifacts repo"
    )
    pull_otifacts_p.set_defaults(which="pull-otifacts")
    # STATUS
    status_p = subp.add_parser(
        "status", help="report the status of a resource (or all resources)"
    )
    status_p.add_argument(
        "resources", nargs="*", help="IDs of the resources to report status on"
    )
    status_p.add_argument(
        "-i", "--ids-only", action="store_true", default=False, help="just list the IDs"
    )
    status_p.add_argument(
        "--by-status",
        action="store_true",
        default=False,
        help="group the report by status",
    )
    status_p.add_argument(
        "--terminal",
        action="store_true",
        default=False,
        help="Report only on the terminalized resource of each type.",
    )
    status_p.set_defaults(which="status")
    # DOWNLOAD
    download_p = subp.add_parser(
        "download", help="download an artifact to your local filesystem"
    )
    download_p.add_argument(
        "resources", nargs="+", help="IDs of the resources to download"
    )
    download_p.set_defaults(which="download")
    # UNPACK
    unpack_p = subp.add_parser(
        "unpack", help="unpack an resource (downloads if necessary)"
    )
    unpack_p.add_argument("resources", nargs="+", help="IDs of the resources to unpack")
    unpack_p.set_defaults(which="unpack")
    # NORMALIZE
    normalize_p = subp.add_parser(
        "normalize", help="converts to the OTT format (unpacks if necessary)"
    )
    normalize_p.add_argument(
        "resources", nargs="+", help="IDs of the resources to normalize"
    )
    normalize_p.set_defaults(which="normalize")
    # PARTITION
    partition_p = subp.add_parser("partition", help="Breaks the resource taxon")
    partition_p.add_argument(
        "--strategy",
        default="hard-coded",
        choices=["hard-coded", "previous"],
        help="Strategy for partitioning the resource. "
        "'hard-coded' relies on ID-mapping stored in the taxalotl "
        "code-base (used for OTT and CoL initial partitions). "
        "The 'previous' strategy uses ID mappings for an external resource "
        "that are found in the source field of OTT.",
    )
    partition_p.add_argument(
        "resources", nargs="+", help="IDs of the resources to partitition"
    )
    _add_level_arg(partition_p)
    partition_p.set_defaults(which="partition")

    # INFO
    info_p = subp.add_parser("info", help="Report statistics about a resource")
    info_p.add_argument("resources", nargs="+", help="IDs of the resources")
    _add_level_arg(info_p)
    info_p.set_defaults(which="info")

    # GREP
    grep_p = subp.add_parser("grep", help="Search the parsed taxonomies")
    grep_p.add_argument("resources", nargs="+", help="IDs of the resources")
    _add_level_arg(grep_p)
    grep_p.add_argument("--name", help="pattern for a name", nargs=1, type=str)
    grep_p.add_argument("--tax-id", help="ID for taxon", nargs=1, type=str)
    grep_p.add_argument(
        "--target",
        default="both",
        choices=["both", "taxa", "synonyms"],
        help="Search in taxa, synonyms or both",
    )
    grep_p.set_defaults(which="grep")

    # ADD-MAPPING
    add_mapping_p = subp.add_parser(
        "add-mapping",
        help="Manually associate an external id to an OTT ID in the partitioned OTT dir",
    )
    add_mapping_p.add_argument(
        "--ott-id", help="The OTT ID to add the mapping", nargs=1, type=str
    )
    add_mapping_p.add_argument(
        "--external-id", help="The OTT ID to add the mapping", nargs=1, type=str
    )
    add_mapping_p.set_defaults(which="add-mapping")

    # CLEAN-PARTITION
    clean_p = subp.add_parser(
        "clean-partition",
        help="remove the results the partition+enforce-new-separator for a resource.",
    )
    clean_p.add_argument("resources", nargs="*", help="IDs of the resources to clean")
    clean_p.set_defaults(which="clean-partition")

    # Handle --show-completions differently from the others, because
    #   argparse does not help us out here... at all
    if "--show-completions" in sys.argv:
        a = sys.argv[1:]
        univ = frozenset(
            [
                "--config",
            ]
        )
        sel_cmd = None
        for c in all_cmds:
            if c in a:
                if sel_cmd is None:
                    sel_cmd = c
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
        elif sel_cmd in res_dep_cmds:
            comp_list = _cmd_completion(a, sel_cmd)

        sys.stdout.write("{}\n".format(" ".join(comp_list)))
    else:
        rc = main_post_parse(p.parse_args())
        sys.exit(rc)


def _cmd_completion(arg_list, sel_cmd):
    a = arg_list

    # From Ned Batchelder's answer on http://stackoverflow.com/a/14728477
    class ArgumentParserError(Exception):
        pass

    # noinspection PyClassHasNoInit
    class ThrowingArgumentParser(argparse.ArgumentParser):
        def error(self, message):
            raise ArgumentParserError(message)

    fake_parser = ThrowingArgumentParser()
    fake_parser.add_argument("--config", type=str)
    fake_parser.add_argument("blah", nargs="*")
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
                taxalotl_config.resources_mgr.abstract_input_resource_types()
            )
    except Exception as _excep:
        _LOG.warning("Exception: {}".format(_excep))
        pass

    if sel_cmd == "status":
        if "-i" not in a and "--ids-only" not in a:
            comp_list.extend(["-i", "--ids-only"])
        for x in ["--by-status", "--terminal"]:
            if x not in a:
                comp_list.extend([x])
    elif sel_cmd == "grep":
        if "--name" == a[-1] or (len(a) > 1 and "--name" == a[-2]):
            comp_list = []
        else:
            acl = ["--level", "--strategy"]
            comp_list = _add_level_and_other_completions(a, comp_list, acl)
    elif sel_cmd == "add-mapping":
        arg_comp_list = ["--ott-id", "--external-id"]
        found = False
        for ac in arg_comp_list:
            if ac == a[-1] or (len(a) > 1 and ac == a[-2]):
                found = True
                comp_list = []
        if not found:
            for x in arg_comp_list:
                if x not in a:
                    comp_list.extend([x])
    elif sel_cmd == "partition":
        if "--strategy" == a[-1] or (len(a) > 1 and "--strategy" == a[-2]):
            comp_list = list(["hard-coded", "previous"])
        else:
            acl = ["--level", "--strategy"]
            comp_list = _add_level_and_other_completions(a, comp_list, acl)
    return comp_list


def _add_level_and_other_completions(arg_list, comp_list, arg_comp_list):
    a = arg_list
    if "--level" == a[-1] or (len(a) > 1 and "--level" == a[-2]):
        comp_list = list(NONTERMINAL_PART_NAMES)
    else:
        for x in arg_comp_list:
            if x not in a:
                comp_list.extend([x])
    return comp_list


if __name__ == "__main__":
    main()
