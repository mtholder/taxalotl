#!/usr/bin/env python
from __future__ import print_function
import sys
import os
from taxalotl import TaxalotlConfig
from peyotl import get_logger
_LOG = get_logger(__name__)

def terminalize_if_needed(wrapper, action):
    if wrapper.is_abstract:
        orig_w = wrapper
        wrapper = wrapper.get_leaf_obj()
        if wrapper is not orig_w:
            m = "{} action being performed on {} as the most recent version of {}"
            _LOG.info(m.format(action, wrapper.id, orig_w.id))
    return wrapper

def download_resource(taxalotl_config, id_list):
    rm = taxalotl_config.resources_mgr
    res = rm.resources
    print(res.keys())
    for rid in id_list:
        try:
            rw = res[rid]
        except:
            raise ValueError("Unknown resource ID '{}'".format(rid))
        rw = terminalize_if_needed(rw, 'download')
        if rw.has_been_downloaded(taxalotl_config):
            m = "{} was already present at {}"
            _LOG.info(m.format(rw.id, rw.download_filepath(taxalotl_config)))
        else:
            rw.download(taxalotl_config)

if __name__ == "__main__":
    import argparse
    description = "The main CLI for taxalotl"
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--resources-dir", type=str, help="the resources directory (optional)")
    p.add_argument("--config", type=str, help="the taxalotl.conf filepath (optional)")
    p.set_defaults(which="all")
    subp = p.add_subparsers(help="command help")
    download_p = subp.add_parser('download', help="download an artifact to your local filesystem")
    download_p.add_argument('resources', nargs="+", help="IDs of the resources to download")
    download_p.set_defaults(which="download")
    unpack_p = subp.add_parser('unpack', help="unpack an artifact to your local filesystem (downloads if necessary)")
    unpack_p.add_argument('resources', nargs="+", help="IDs of the resources to download")
    unpack_p.set_defaults(which="unpack")
    args = p.parse_args()
    taxalotl_config = TaxalotlConfig(filepath=args.config, resources_dir=args.resources_dir)
    if args.which == 'download':
        download_resource(taxalotl_config, args.resources)