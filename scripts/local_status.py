#!/usr/bin/env python
from __future__ import print_function
import sys
from taxalotl import TaxalotlConfig
def main(taxalotl_config):
    rm = taxalotl_config.resources_mgr
    ks = rm.resources.keys()
    ks.sort()
    for k in ks:
        r = rm.resources[k]
        r.write_status(sys.stdout, taxalotl_config, indent='  ')


if __name__ == "__main__":
    import argparse
    description = "Reports on the known resources and whether or not you have downloaded an processed them."
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--resources-dir", type=str, help="the resources directory (optional)")
    p.add_argument("--config", type=str, help="the taxalotl.conf filepath (optional)")
    args = p.parse_args()
    tc = TaxalotlConfig(filepath=args.config, resources_dir=args.resources_dir)
    main(tc)