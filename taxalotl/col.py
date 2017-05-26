from __future__ import print_function

import codecs

from peyotl import (get_logger, shorter_fp_form)

from taxalotl.resource_wrapper import ExternalTaxonomyWrapper

_LOG = get_logger(__name__)

COL_PARTMAP = {
    'Archaea': frozenset([33524792]),
    'Bacteria': frozenset([33521420]),
    'Eukaryota': frozenset([33521288, 33521351, 33521293, 33521595, 33523363]),
    'Archaeplastida': frozenset([33521293]),
    'Glaucophyta': frozenset([33531372]),
    'Rhodophyta': frozenset([33531376]),
    'Chloroplastida': frozenset([33521294, 33522836, 33522993, 33531380, 33531384, 33534508]),
    'Fungi': frozenset([33521351]),
    'Metazoa': frozenset([33521288]),
    'Annelida': frozenset([33521477]),
    'Arthropoda': frozenset([33521342]),
    'Malacostraca': frozenset([33521356]),
    'Arachnida': frozenset([33521366]),
    'Insecta': frozenset([33521474]),
    'Diptera': frozenset([33521519]),
    'Coleoptera': frozenset([33521475]),
    'Lepidoptera': frozenset([33521695]),
    'Hymenoptera': frozenset([33521644]),
    'Bryozoa': frozenset([33524015]),
    'Chordata': frozenset([33521289]),
    'Cnidaria': frozenset([33522061]),
    'Ctenophora': frozenset([33521313]),
    'Mollusca': frozenset([33521301]),
    'Nematoda': frozenset([33526516]),
    'Platyhelminthes': frozenset([33521309]),
    'Porifera': frozenset([33527549]),
    'Viruses': frozenset([33521407]),
}


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
    complete_taxon_fp = tax_part.tax_fp
    syn_fp = tax_part.syn_fp
    assert not syn_fp
    syn_by_id = tax_part.syn_by_id
    ptp = shorter_fp_form(complete_taxon_fp)
    with codecs.open(complete_taxon_fp, 'rU', encoding='utf-8') as inp:
        iinp = iter(inp)
        tax_part.taxon_header = iinp.next()
        prev_line = None
        # vt = unicode('\x0b') # Do some lines have vertical tabs? Of course they do....
        # istwo = unicode('\x1e')
        for n, line in enumerate(iinp):
            if not line.endswith('\n'):
                if prev_line:
                    prev_line = prev_line + line[:-1]
                else:
                    prev_line = line[:-1]
                continue
            elif prev_line:
                line = prev_line + line
                prev_line = ''
            ls = line.split('\t')
            if n % 1000 == 0:
                _LOG.info(' read taxon {} from {}'.format(n, ptp))
            try:
                col_id, accept_id, par_id = ls[0], ls[4], ls[5]
                col_id = int(col_id)
                if accept_id:
                    accept_id = int(accept_id)
                    syn_by_id.setdefault(accept_id, []).append((col_id, line))
                else:
                    tax_part.read_taxon_line(col_id, par_id, line)
            except:
                _LOG.exception("Exception parsing line {}:\n{}".format(1 + n, line))
                raise


# noinspection PyAbstractClass
class CoLExternalTaxonomyWrapper(ExternalTaxonomyWrapper):
    taxon_filename = 'taxa.txt'
    synonyms_filename = None
    partition_parsing_fn = staticmethod(partition_col_by_root_id)

    def __init__(self, obj, parent=None, refs=None):
        ExternalTaxonomyWrapper.__init__(self, obj, parent=parent, refs=refs)

    @property
    def partition_source_dir(self):
        return self.unpacked_filepath

    def get_primary_partition_map(self):
        return COL_PARTMAP
