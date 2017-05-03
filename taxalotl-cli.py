#!/usr/bin/env python
from __future__ import print_function
import sys
from taxalotl import TaxalotlConfig
from peyotl import get_logger

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
    for i in id_list:
        r = res[i]
        if not r.has_been_downloaded():
            if r.is_abstract:
                a_list.append(i)
            else:
                nd_list.append(i)
        elif not r.has_been_unpacked():
            dnu_list.append(i)
        elif not r.has_been_normalized():
            unn_list.append(i)
        else:
            n_list.append(i)
    return [["abstract classes", a_list],
            ["not downloaded", nd_list],
            ["downloaded, but not unpacked", dnu_list],
            ["unpacked, but not normalized", unn_list],
            ["normalized", n_list],
           ]

def status_of_resources(taxalotl_config,
                        id_list,
                        ids_only=False,
                        by_status=False):
    terminalize = True
    res = taxalotl_config.resources_mgr.resources
    if not id_list:
        terminalize = False
        id_list = list(res.keys())
        id_list.sort()
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
                rw = taxalotl_config.get_terminalized_res_by_id(rid, 'status')
            else:
                rw = taxalotl_config.get_resource_by_id(rid)
            rw.write_status(out_stream, indent=' '*4)


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


def main(args):
    taxalotl_config = TaxalotlConfig(filepath=args.config, resources_dir=args.resources_dir)
    try:
        if args.which == 'download':
            download_resources(taxalotl_config, args.resources)
        elif args.which == 'status':
            status_of_resources(taxalotl_config,
                                args.resources,
                                ids_only=args.ids_only,
                                by_status=args.by_status)
        elif args.which == 'unpack':
            unpack_resources(taxalotl_config, args.resources)
        elif args.which == 'normalize':
            normalize_resources(taxalotl_config, args.resources)
    except Exception as x:
        if taxalotl_config.crash_with_stacktraces:
            raise
        sys.exit('taxalotl-cli: Exiting with exception:\n{}'.format(x))


if __name__ == "__main__":
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
    subp = p.add_subparsers( help="command help")
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
    normalize_p.add_argument('resources', nargs="+", help="IDs of the resources to unpack")
    normalize_p.set_defaults(which="normalize")
    # Handle --show-completions differently from the others, because
    #   argparse does not help us out here... at all
    if "--show-completions" in sys.argv:
        a = sys.argv[1:]
        univ = frozenset(['--resources-dir', '--config'])
        res_dep_cmds = ['status', 'download', 'unpack', 'normalize']
        all_cmds = res_dep_cmds
        sel_cmd = None
        for c in all_cmds:
            if c in a:
                sel_cmd = c # perhaps should worry about multiple commands?
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
                class ThrowingArgumentParser(argparse.ArgumentParser):
                    def error(self, message):
                        raise ArgumentParserError(message)
                fake_parser = ThrowingArgumentParser()
                fake_parser.add_argument("--resources-dir", type=str)
                fake_parser.add_argument("--config", type=str)
                fake_parser.add_argument('blah', nargs="*")
                resdir, config = None, None
                comp_list = []
                try:
                    fa, unk = fake_parser.parse_known_args()
                    resdir, config = fa.resources_dir, fa.config
                    taxalotl_config = TaxalotlConfig(filepath=config, resources_dir=resdir)
                    comp_list = list(taxalotl_config.resources_mgr.resources.keys())
                except:
                    raise
                if sel_cmd == 'status':
                    if '-i' not in a and '--ids-only' not in a:
                        comp_list.extend(["-i","--ids-only"])
                    if '--by-status' not in a:
                        comp_list.extend('--by-status')
        sys.stdout.write('{}\n'.format(' '.join(comp_list)))
    else:
        main(p.parse_args())
