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

def get_res_by_id(resource_dict, res_id):
    try:
        return resource_dict[res_id]
    except:
        raise ValueError("Unknown resource ID '{}'".format(res_id))

def get_terminalized_res_by_id(resource_dict, res_id, for_action):
    rw = get_res_by_id(resource_dict=resource_dict, res_id=res_id)
    return terminalize_if_needed(rw, for_action)


def download_resources(taxalotl_config, id_list):
    res = taxalotl_config.resources_mgr.resources
    for rid in id_list:
        rw = get_terminalized_res_by_id(res, rid, 'download')
        if rw.has_been_downloaded(taxalotl_config):
            m = "{} was already present at {}"
            _LOG.info(m.format(rw.id, rw.download_filepath(taxalotl_config)))
        else:
            rw.download(taxalotl_config)

def status_of_resources(taxalotl_config, id_list):
    res = taxalotl_config.resources_mgr.resources
    terminalize = True
    if not id_list:
        terminalize = False
        id_list = list(res.keys())
        id_list.sort()
    for rid in id_list:
        if terminalize:
            rw = get_terminalized_res_by_id(res, rid, 'status')
        else:
            rw = get_res_by_id(res, rid)
        rw.write_status(sys.stdout, taxalotl_config, indent='  ')


def unpack_resources(taxalotl_config, id_list):
    res = taxalotl_config.resources_mgr.resources
    for rid in id_list:
        rw = get_terminalized_res_by_id(res, rid, 'unpack')
        if not rw.has_been_downloaded(taxalotl_config):
            m = "{} will be downloaded first..."
            _LOG.info(m.format(rw.id))
            download_resource(taxalotl_config, [rw.id])
        if rw.has_been_unpacked(taxalotl_config):
            m = "{} was already present at {}"
            _LOG.info(m.format(rw.id, rw.unpacked_filepath(taxalotl_config)))
        else:
            rw.unpack(taxalotl_config)


if __name__ == "__main__":
    import argparse
    description = "The main CLI for taxalotl"
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--resources-dir", type=str, help="the resources directory (optional)")
    p.add_argument("--config", type=str, help="the taxalotl.conf filepath (optional)")
    p.set_defaults(which="all")
    subp = p.add_subparsers(help="command help")
    #
    download_p = subp.add_parser('download', help="download an artifact to your local filesystem")
    download_p.add_argument('resources', nargs="+", help="IDs of the resources to download")
    download_p.set_defaults(which="download")
    #
    status_p = subp.add_parser('status',
                               help="report the status of a resource (or all resources)")
    status_p.add_argument('resources', nargs="*", help="IDs of the resources to report status on")
    status_p.set_defaults(which="status")
    #
    unpack_p = subp.add_parser('unpack',
                               help="unpack an artifact to your local filesystem (downloads if necessary)")
    unpack_p.add_argument('resources', nargs="+", help="IDs of the resources to unpack")
    unpack_p.set_defaults(which="unpack")

    args = p.parse_args()
    taxalotl_config = TaxalotlConfig(filepath=args.config, resources_dir=args.resources_dir)
    if args.which == 'download':
        download_resources(taxalotl_config, args.resources)
    elif args.which == 'status':
        status_of_resources(taxalotl_config, args.resources)
    elif args.which == 'unpack':
        unpack_resources(taxalotl_config, args.resources)