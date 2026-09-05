#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
import io
import os
import re

# noinspection PyPep8Naming
import xml.etree.ElementTree as ET
import logging

from peyutil import assure_dir_exists

from ..ott_schema import InterimTaxonomyData
from ..resource_wrapper import TaxonomyWrapper
from ..util import OutFile

_LOG = logging.getLogger(__name__)


def normalize_coldp_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    sys.exit("early")
    manifest_fp = os.path.join(source, "meta.xml")
    manifest_root = ET.parse(manifest_fp).getroot()
    core_paths = []
    field2index = {}
    for el in manifest_root.findall("{http://rs.tdwg.org/dwc/text/}core"):
        for sub in el:
            if sub.tag.endswith("}id"):
                field2index["id"] = int(sub.attrib["index"])
            elif sub.tag.endswith("}field"):
                nns = os.path.split(sub.attrib["term"])[-1]
                field2index[nns] = int(sub.attrib["index"])
        for f in el.findall("{http://rs.tdwg.org/dwc/text/}files"):
            for loc in f.findall("{http://rs.tdwg.org/dwc/text/}location"):
                core_paths.append(loc.text.strip())
    if len(core_paths) != 1:
        raise ValueError(
            'Did not find a single core path in DwC file ("{}") found: {}'.format(
                manifest_fp, core_paths
            )
        )
    taxon_fn = core_paths[0]
    proj_out = os.path.join(destination, "projection.tsv")
    if not os.path.exists(proj_out):
        proj_in = os.path.join(source, taxon_fn)
        write_gbif_projection_file(proj_in, proj_out, field2index)
    homemade = {
        "id": 0,
        "parentNameUsageID": 1,
        "acceptedNameUsageID": 2,
        "canonicalName": 3,
        "taxonRank": 4,
        "taxonomicStatus": 5,
        "nameAccordingTo": 6,
    }

    itd = InterimTaxonomyData()
    to_remove, to_ignore, paleos = read_gbif_projection(
        proj_out, itd, homemade, do_gbif_checks=isinstance(res_wrapper, GBIFWrapper)
    )
    add_fake_root(itd)
    remove_if_tips(itd, to_remove)
    o_to_ignore = find_orphaned(itd)
    to_ignore.update(o_to_ignore)
    prune_ignored(itd, to_ignore)
    _LOG.info("writing {} paleodb ids".format(len(paleos)))
    with OutFile(os.path.join(destination, "paleo.tsv")) as paleofile:
        for taxon_id in paleos:
            paleofile.write("{}\n".format(taxon_id))
    res_wrapper.post_process_interim_tax_data(itd)
    itd.write_to_dir(destination)
