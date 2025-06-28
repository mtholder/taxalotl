from __future__ import print_function

import io
import logging

from peyutil import shorter_fp_form

from ..resource_wrapper import TaxonomyWrapper
from ..parsing.darwin_core import normalize_darwin_core_taxonomy

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


# noinspection PyUnreachableCode
def partition_col_by_root_id(tax_part):  # type (TaxonPartition) -> None
    """Reads the serialized taxonomy of the parent, adds the easy lines to their partition element,
    and returns dicts needed to finish the assignments.

    Signature for partition functions. Takes:
      1. abs path of taxonomy file for parent taxon
      2. list of PartitionElements whose roots are sets that specify IDs that are the
        roots of the subtrees that are to go in each partition elemen.
    Returns a tuple:
      0. par_id ->[child_id] dict,
      1. id -> partition_element dict for already assigned IDs,
      2. id -> line dict - may only have unassigned IDs in it,
      3. synonym id -> [(accepted_id, line), ] for any synonyms
      4. roots_set - a frozen set of the union of the partition element roots
      5. the rootless partition element ("garbage_bin" for all unassigned IDs)
      6. header for taxon file
      7. header for synonyms file (or None)

    """
    assert False
    complete_taxon_fp = tax_part.tax_fp
    syn_fp = tax_part.input_synonyms_filepath
    assert not syn_fp
    syn_by_id = tax_part._syn_by_id
    ptp = shorter_fp_form(complete_taxon_fp)
    with io.open(complete_taxon_fp, "r", encoding="utf-8") as inp:
        iinp = iter(inp)
        tax_part.taxon_header = next(iinp)
        prev_line = None
        # vt = unicode('\x0b') # Do some lines have vertical tabs? Of course they do....
        # istwo = unicode('\x1e')
        for n, line in enumerate(iinp):
            if not line.endswith("\n"):
                if prev_line:
                    prev_line = prev_line + line[:-1]
                else:
                    prev_line = line[:-1]
                continue
            elif prev_line:
                line = prev_line + line
                prev_line = ""
            ls = line.split("\t")
            if n % 1000 == 0:
                _LOG.info(" read taxon {} from {}".format(n, ptp))
            try:
                col_id, accept_id, par_id = ls[0], ls[4], ls[5]
                col_id = int(col_id)
                if accept_id:
                    try:
                        accept_id = int(accept_id)
                    except:
                        if n == 0:
                            continue
                    syn_by_id.setdefault(accept_id, []).append((col_id, line))
                else:
                    tax_part.read_taxon_line(col_id, par_id, line)
            except Exception:
                _LOG.exception("Exception parsing line {}:\n{}".format(1 + n, line))
                raise


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

    def _post_process_tree(self, tree):
        self.collapse_incertae_sedis_by_name_prefix(tree, "not assigned")

    def post_process_interim_tax_data(self, interim_tax_data):
        self.collapse_as_incertae_sedis_interim_tax_data(
            interim_tax_data, "not assigned"
        )
