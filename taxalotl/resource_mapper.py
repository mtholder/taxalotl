#!/usr/bin/env python
"""Should not be imported directly! (to avoid circular dependencies)
Imported by resource manager funtions.

Holds mapping of base_id to class of wrapper to instantiate.
"""
from taxalotl.resource_wrapper import ExternalTaxonomyWrapper
from taxalotl.ott import OTTaxonomyWrapper, OTTaxonomyIdListWrapper
from taxalotl.silva import SilvaIdListWrapper, SilvaWrapper
from taxalotl.col import CoLExternalTaxonomyWrapper
from taxalotl.irmng import IRMNGWrapper
from taxalotl.darwin_core import GBIFWrapper
from taxalotl.ncbi import NCBIWrapper

BASE_ID_TO_RES_TYPE = {
    'col': CoLExternalTaxonomyWrapper,
    'gbif': GBIFWrapper,
    'irmng': IRMNGWrapper,
    'ncbi': NCBIWrapper,
    'ott': OTTaxonomyWrapper,
    'silva': SilvaWrapper,
}

wrapper_types = set()


# noinspection PyTypeChecker
def _fill_wrapper():
    global wrapper_types
    wrapper_types.update(BASE_ID_TO_RES_TYPE.values())
    wrapper_types.add(ExternalTaxonomyWrapper)
    wrapper_types.add(SilvaIdListWrapper)
    wrapper_types.add(OTTaxonomyIdListWrapper)


_fill_wrapper()