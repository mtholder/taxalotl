#!/usr/bin/env python3

import os
import xml.etree.ElementTree as ElementTree
from peyotl import get_logger

_LOG = get_logger(__name__)


def _get_col_fields_and_taxa_fp(col_dir):
    manifest_fp = os.path.join(col_dir, 'meta.xml')
    manifest_root = ElementTree.parse(manifest_fp).getroot()
    core_paths = []
    field2index = {}
    for el in manifest_root.findall('{http://rs.tdwg.org/dwc/text/}core'):
        for sub in el:
            if sub.tag.endswith('}id'):
                field2index['id'] = int(sub.attrib['index'])
            elif sub.tag.endswith('}field'):
                nns = os.path.split(sub.attrib['term'])[-1]
                field2index[nns] = int(sub.attrib['index'])
        for f in el.findall('{http://rs.tdwg.org/dwc/text/}files'):
            for loc in f.findall('{http://rs.tdwg.org/dwc/text/}location'):
                core_paths.append(loc.text.strip())
    if len(core_paths) != 1:
        raise ValueError(
            'Did not find a single core path in DwC file ("{}") found: {}'.format(
                manifest_fp, core_paths))
    taxon_fn = core_paths[0]
    return field2index, os.path.join(col_dir, taxon_fn)


class CoLDumpParser(object):
    def __init__(self, col_dir):
        f2i, taxon_fp = _get_col_fields_and_taxa_fp(col_dir)
        ts = [(v, k) for k, v in f2i.items() if k != 'id']
        ts.sort()
        self._col_order = [i[1] for i in ts]
        self._taxon_fp = taxon_fp
        self.taxonID_idx = f2i["taxonID"]
        self.identifier_idx = f2i["identifier"]
        self.datasetID_idx = f2i["datasetID"]
        self.datasetName_idx = f2i["datasetName"]
        self.acceptedNameUsageID_idx = f2i["acceptedNameUsageID"]
        self.parentNameUsageID_idx = f2i["parentNameUsageID"]
        self.taxonomicStatus_idx = f2i["taxonomicStatus"]
        self.taxonRank_idx = f2i["taxonRank"]
        self.verbatimTaxonRank_idx = f2i["verbatimTaxonRank"]
        self.scientificName_idx = f2i["scientificName"]
        self.kingdom_idx = f2i["kingdom"]
        self.phylum_idx = f2i["phylum"]
        self.class_idx = f2i["class"]
        self.order_idx = f2i["order"]
        self.superfamily_idx = f2i["superfamily"]
        self.family_idx = f2i["family"]
        self.genericName_idx = f2i["genericName"]
        self.genus_idx = f2i["genus"]
        self.subgenus_idx = f2i["subgenus"]
        self.specificEpithet_idx = f2i["specificEpithet"]
        self.infraspecificEpithet_idx = f2i["infraspecificEpithet"]
        self.scientificNameAuthorship_idx = f2i["scientificNameAuthorship"]
        self.source_idx = f2i["source"]
        self.namePublishedIn_idx = f2i["namePublishedIn"]
        self.nameAccordingTo_idx = f2i["nameAccordingTo"]
        self.modified_idx = f2i["modified"]
        self.description_idx = f2i["description"]
        self.taxonConceptID_idx = f2i["taxonConceptID"]
        self.scientificNameID_idx = f2i["scientificNameID"]
        self.references_idx = f2i["references"]
        self.isExtinct_idx = f2i["isExtinct"]

    def gen_record_lists(self):
        """Checks header then is a generator for records split as lists of strings"""
        with open(self._taxon_fp, 'rU', encoding='utf-8') as inp:
            line_iter = iter(inp)
            headers = next(line_iter)[:-1].split('\t')
            if headers[1:] != self._col_order[1:]:
                m = 'Expecting headers to be {}, but found {}'
                raise ValueError(m.format(self._col_order, headers))
            for n, line in enumerate(inp):
                yield line[:-1].split('\t')

    def gen_records(self):
        """Generator for records, processes int and bool entries into Python types"""
        int_attribs = [self.taxonID_idx,
                       self.acceptedNameUsageID_idx,
                       self.parentNameUsageID_idx]
        bool_attribs = [self.isExtinct_idx]
        typed_attribs = int_attribs + bool_attribs
        str_attribs = [i for i in range(len(self._col_order)) if i not in typed_attribs]
        type_err_template = 'Expecting {} in column index {} ({}) but found "{}" in line {} ("{}")'
        for n, sl in enumerate(self.gen_record_lists()):
            try:
                for i in int_attribs:
                    v = sl[i]
                    if not v:
                        sl[i] = None
                    else:
                        try:
                            sl[i] = int(v)
                        except Exception:
                            m = type_err_template.format(i, self._col_order[i], v, 2 + n, sl)
                            raise RuntimeError(m)
                for i in bool_attribs:
                    v = sl[i].lower()
                    if not v:
                        sl[i] = None
                    else:
                        try:
                            assert v == 'true' or v == 'false'
                            sl[i] = v == 'true'
                        except Exception:
                            m = type_err_template.format(i, self._col_order[i], v, 2 + n, sl)
                            raise RuntimeError(m)
                for i in str_attribs:
                    if not sl[i]:
                        sl[i] = None
                # if sl[self.taxonID_idx] == 38525922:
                #    pass
                yield sl
            except Exception:
                _LOG.warn('Exception parsing line:\n{}'.format(sl))
                raise


class ParseCondenser(object):

    def __init__(self, col_dump_parser):
        p = col_dump_parser
        self.rec_gen = p
        self._taxon_fp = p._taxon_fp
        retained_fields = [(p.acceptedNameUsageID_idx, "acceptedNameUsageID_idx"),
                           (p.datasetID_idx, "datasetID_idx"),
                           (p.description_idx, "description_idx"),
                           (p.genericName_idx, "genericName_idx"),
                           (p.identifier_idx, "identifier_idx"),
                           (p.infraspecificEpithet_idx, "infraspecificEpithet_idx"),
                           (p.isExtinct_idx, "isExtinct_idx"),
                           (p.modified_idx, "modified_idx"),
                           (p.nameAccordingTo_idx, "nameAccordingTo_idx"),
                           (p.parentNameUsageID_idx, "parentNameUsageID_idx"),
                           (p.references_idx, "references_idx"),
                           (p.scientificName_idx, "scientificName_idx"),
                           (p.scientificNameAuthorship_idx, "scientificNameAuthorship_idx"),
                           (p.scientificNameID_idx, "scientificNameID_idx"),
                           (p.specificEpithet_idx, "specificEpithet_idx"),
                           (p.subgenus_idx, "subgenus_idx"),
                           (p.taxonConceptID_idx, "taxonConceptID_idx"),
                           (p.taxonID_idx, "taxonID_idx"),
                           (p.taxonomicStatus_idx, "taxonomicStatus_idx"),
                           (p.taxonRank_idx, "taxonRank_idx"),
                           (p.verbatimTaxonRank_idx, "verbatimTaxonRank_idx"),
                           ]
        retained_fields.sort()
        self._col_order = []
        self._source_idx = []
        for condensed_idx, t in enumerate(retained_fields):
            raw_idx, attr_name = t
            assert attr_name.endswith('_idx')
            col_name = attr_name[:-4]
            self._col_order.append(col_name)
            self._source_idx.append(raw_idx)
            setattr(self, attr_name, condensed_idx)
        transformed_index = len(self._source_idx)
        self.transformed_columns = ('rank_sort_num_idx',  # integer to make unranked taxa sortable..
                                    'taxon_supp_inf_idx',  # slot for taxon_supp_inf key
                                    'entry_trav_int_idx',  # lowest traversal # in clade
                                    'exit_trav_int_idx',  # highest traversal # in clade
                                    )
        for offset, attr_name in enumerate(self.transformed_columns):
            setattr(self, attr_name, transformed_index + offset)
        self._extend_arg = [None] * len(self.transformed_columns)

    def gen_records(self):
        for raw in self.rec_gen.gen_records():
            condensed = [raw[i] for i in self._source_idx]
            condensed.extend(self._extend_arg)
            yield condensed
