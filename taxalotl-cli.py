#!/usr/bin/env python
from __future__ import print_function
import sys
from taxalotl import TaxalotlConfig
from peyotl import get_logger

_LOG = get_logger(__name__)


def download_resources(taxalotl_config, id_list):
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, 'download')
        if rw.has_been_downloaded():
            m = "{} was already present at {}"
            _LOG.info(m.format(rw.id, rw.download_filepath(taxalotl_config)))
        else:
            rw.download(taxalotl_config)


def status_of_resources(taxalotl_config, id_list):
    terminalize = True
    if not id_list:
        terminalize = False
        id_list = list(taxalotl_config.resources_mgr.resources.keys())
        id_list.sort()
    for rid in id_list:
        if terminalize:
            rw = taxalotl_config.get_terminalized_res_by_id(rid, 'status')
        else:
            rw = taxalotl_config.get_resource_by_id(rid)
        rw.write_status(sys.stdout, indent='  ')


def unpack_resources(taxalotl_config, id_list):
    for rid in id_list:
        rw = taxalotl_config.get_terminalized_res_by_id(rid, 'unpack')
        if not rw.has_been_downloaded():
            m = "{} will be downloaded first..."
            _LOG.info(m.format(rw.id))
            download_resources(taxalotl_config, [rw.id])
        if rw.has_been_unpacked():
            m = "{} was already present at {}"
            _LOG.info(m.format(rw.id, rw.unpacked_filepath()))
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
            _LOG.info(m.format(rw.id, rw.normalized_filepath()))
        else:
            rw.normalize()


def main(args):
    taxalotl_config = TaxalotlConfig(filepath=args.config, resources_dir=args.resources_dir)
    if args.which == 'download':
        download_resources(taxalotl_config, args.resources)
    elif args.which == 'status':
        status_of_resources(taxalotl_config, args.resources)
    elif args.which == 'unpack':
        unpack_resources(taxalotl_config, args.resources)
    elif args.which == 'normalize':
        normalize_resources(taxalotl_config, args.resources)


if __name__ == "__main__":
    import argparse

    description = "The main CLI for taxalotl"
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--resources-dir", type=str, help="the resources directory (optional)")
    p.add_argument("--config", type=str, help="the taxalotl.conf filepath (optional)")
    p.set_defaults(which="all")
    subp = p.add_subparsers(help="command help")
    # STATUS
    status_p = subp.add_parser('status',
                               help="report the status of a resource (or all resources)")
    status_p.add_argument('resources', nargs="*", help="IDs of the resources to report status on")
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
    # call main
    main(p.parse_args())
