#!/usr/bin/env python
"""Should not be imported directly! (to avoid circular dependencies)
Imported by resource manager funtions.

Holds mapping of base_id to class of wrapper to instantiate.
"""
from taxalotl.parsing.col import CoLTaxonomyWrapper
from taxalotl.parsing.darwin_core import GBIFWrapper
from taxalotl.parsing.irmng import IRMNGWrapper
from taxalotl.parsing.ncbi import NCBIWrapper
from taxalotl.parsing.ott import OTTaxonomyWrapper, OTTaxonomyIdListWrapper
from taxalotl.parsing.plant_list import PlantListWrapper
from taxalotl.parsing.silva import SilvaIdListWrapper, SilvaWrapper, SilvaToNCBIMappingListWrapper
from taxalotl.resource_wrapper import TaxonomyWrapper
from taxalotl.parsing.worms import WormsTaxonomyWrapper

BASE_ID_TO_RES_TYPE = {
    'col': CoLTaxonomyWrapper,
    'gbif': GBIFWrapper,
    'irmng_raw': IRMNGWrapper,
    'irmng_ot': TaxonomyWrapper,
    'ncbi': NCBIWrapper,
    'ott': OTTaxonomyWrapper,
    'silva': SilvaWrapper,
    'plantlist': PlantListWrapper,
    'worms': WormsTaxonomyWrapper,
}

wrapper_types = set()


# noinspection PyTypeChecker
def _fill_wrapper():
    global wrapper_types
    wrapper_types.update(BASE_ID_TO_RES_TYPE.values())
    wrapper_types.add(TaxonomyWrapper)
    wrapper_types.add(SilvaIdListWrapper)
    wrapper_types.add(SilvaToNCBIMappingListWrapper)
    wrapper_types.add(OTTaxonomyIdListWrapper)


_fill_wrapper()
