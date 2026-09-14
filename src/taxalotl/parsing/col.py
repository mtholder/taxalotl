from __future__ import print_function

import io
import logging

from peyutil import shorter_fp_form

from ..resource_wrapper import TaxonomyWrapper
from ..parsing.darwin_core import (
    GBIFWrapper,
    normalize_darwin_core_taxonomy,
)
from ..parsing.coldp import normalize_coldp_taxonomy

_LOG = logging.getLogger(__name__)

COL_PARTMAP = {
    "Archaea": frozenset([52435722]),
    "Bacteria": frozenset([52433432]),
    "Eukaryota": frozenset([52433499, 52435027, 52433974, 52433370]),
    "Archaeplastida": frozenset([52433499]),
    "Fungi": frozenset([52433393]),
    "Metazoa": frozenset([52433370]),
    "Viruses": frozenset([52433426]),
    "Glaucophyta": frozenset([52444130]),
    "Rhodophyta": frozenset([52444134]),
    "Chloroplastida": frozenset(
        [
            52442327,
            52442210,
            52442148,
            52434330,
            52434201,
            52433500,
        ]
    ),
    "Annelida": frozenset([52433489]),
    "Arthropoda": frozenset([52433375]),
    "Malacostraca": frozenset([52433389]),
    "Arachnida": frozenset([52433402]),
    "Insecta": frozenset([52433376]),
    "Diptera": frozenset([52433521]),
    "Coleoptera": frozenset([52433486]),
    "Lepidoptera": frozenset([52433663]),
    "Hymenoptera": frozenset([52433621]),
    "Bryozoa": frozenset([52442814]),
    "Chordata": frozenset([52433371]),
    "Cnidaria": frozenset([52433398]),
    "Ctenophora": frozenset([52443092]),
    "Mollusca": frozenset([52440786]),
    "Nematoda": frozenset([52436787]),
    "Platyhelminthes": frozenset([52443117]),
    "Porifera": frozenset([52442836]),
}


# noinspection PyAbstractClass
class CoLTaxonomyWrapper(TaxonomyWrapper):
    taxon_filename = "taxonomy.tsv"
    # synonyms_filename = None
    # partition_parsing_fn = staticmethod(partition_col_by_root_id)
    schema = {"http://rs.tdwg.org/dwc/"}

    def __init__(self, obj, parent=None, refs=None):
        TaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    @property
    def partition_source_dir(self):
        return self.normalized_filedir

    def get_primary_partition_map(self):
        return COL_PARTMAP

    def normalize(self):
        normalize_darwin_core_taxonomy(
            self.unpacked_filepath, self.normalized_filedir, self
        )

    def post_process_interim_tax_data(self, interim_tax_data):
        self.collapse_as_incertae_sedis_interim_tax_data(
            interim_tax_data, "not assigned"
        )


class CoLXRTaxonomyWrapper(TaxonomyWrapper):
    taxon_filename = "NameUsage.tsv"
    schema = {"https://github.com/CatalogueOfLife/coldp/releases/tag/v1.2.0"}

    def normalize(self):
        normalize_coldp_taxonomy(self.unpacked_filepath, self.normalized_filedir, self)
