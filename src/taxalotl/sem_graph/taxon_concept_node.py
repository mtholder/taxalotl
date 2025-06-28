#!/usr/bin/env python
from .graph_node import SemGraphNode

KNOWN_FLAGS = frozenset(["hidden", "sibling_higher", "extinct"])


class TaxonConceptSemNode(SemGraphNode):
    def __init__(self, sem_graph, foreign_id):
        d = {"class_tag": "tc", "context_id": foreign_id}
        super(TaxonConceptSemNode, self).__init__(sem_graph, d)
        self.is_child_of = None
        self.rank = None
        self.has_name = None
        self.undescribed = None
        self._is_synonym_of = None
        self.problematic_synonyms = None
        self.synonyms = None
        self.syn_type = None
        self.former_ranks = None
        self.hybrid = None
        self.incertae_sedis = None
        self.other_flags = None
        self._child_set = None
        self.mapped_to = None

    def add_verbatim_name(self, name_str):
        n = self.graph._add_verbatim_name(self, name_str)
        self.claim_name(n)
        return n

    @property
    def most_terminal_name(self):
        if self.has_name is None:
            return None
        return self.has_name.most_terminal_name

    @property
    def problematic_synonym_list(self):
        return self.problematic_synonyms if self.problematic_synonyms else []

    @property
    def synonym_list(self):
        return self.synonyms if self.synonyms else []

    @property
    def child_set(self):
        if self._child_set is None:
            self._child_set = set()
        return self._child_set

    def find_in_subtree(self, test):
        r = []
        if test(self):
            r.append(self)
        for c in self.child_set:
            r.extend(c.find_in_subtree(test))
        return r

    def find_valid_species(self, genus_name, species_name):
        def test_fn(tc):
            return (
                (not tc.undescribed)
                and (tc.rank == "species")
                and tc.is_valid_for_sp_epithet(species_name)
            )

        return self.find_in_subtree(test=test_fn)

    def find_undescribed_species_for_name(self, genus_name, species_name):
        def test_fn(tc):
            return (
                tc.undescribed
                and (tc.rank == "species")
                and tc.is_valid_for_sp_epithet(species_name)
            )

        return self.find_in_subtree(test=test_fn)

    def explain(self, out):
        cn = self.canonical_name.name if self.canonical_name else ""
        out.write('"{}" Taxon({}) '.format(cn, self.canonical_id))
        if self._is_synonym_of:
            out.write('synonym of "{}"'.format(self.valid_name.name))
        else:
            out.write(" rank={}".format(self.rank if self.rank else "?"))

    def is_valid_for_sp_epithet(self, species_name):
        if self._is_synonym_of:
            return False
        if self.has_name.sp_epithet and self.has_name.sp_epithet.name == species_name:
            return True
        for syn in self.synonym_list:
            cn = syn.has_name.sp_epithet
            if cn.name == species_name:
                return True
        return False

    def is_valid_for_name(self, genus_name):
        if self._is_synonym_of:
            return False
        if self.canonical_name.name == genus_name:
            return True
        for syn in self.synonym_list:
            cn = syn.canonical_name
            if cn.name == genus_name:
                return True
        return False

    @property
    def valid_name(self):
        if self._is_synonym_of:
            return self._is_synonym_of.valid_name
        return self.canonical_name

    @property
    def canonical_name(self):
        if not self.has_name:
            return None
        if self.has_name.normalized:
            return self.has_name.normalized
        return self.has_name.canonical_name

    @property
    def is_synonym_of(self):
        return bool(self._is_synonym_of)

    @property
    def is_the_valid_name(self):
        return not self.is_synonym_of

    @property
    def is_specimen_based(self):
        return self.rank in ["species", "subspecies", "infraspecies"]

    def claim_is_child_of(self, par_sem_node):
        assert self.is_child_of is None
        self.is_child_of = par_sem_node
        par_sem_node.child_set.add(self)

    def claim_rank(self, rank):
        assert self.rank is None
        self.rank = rank

    def claim_name(self, name_sem_node):
        assert self.has_name is None
        self.has_name = name_sem_node

    def claim_undescribed(self):
        self.undescribed = True

    def claim_hybrid(self):
        self.hybrid = True

    def claim_incertae_sedis(self):
        self.incertae_sedis = True

    def claim_flag(self, flag):
        if flag == "incertae_sedis":
            self.claim_incertae_sedis()
        elif flag == "hybrid":
            self.claim_hybrid()
        else:
            if flag not in KNOWN_FLAGS:
                raise ValueError('Unknown flag "{}"'.format(flag))
            self._add_other_flag(flag)

    def _add_other_flag(self, flag):
        if self.other_flags:
            if flag in self.other_flags:
                return
            self.other_flags.append(flag)
            self.other_flags.sort()
        else:
            self.other_flags = [flag]

    @property
    def predicates(self):
        return [
            "hybrid",
            "is_child_of",
            "rank",
            "has_name",
            "id",
            "undescribed",
            "is_synonym",
            "incertae_sedis",
            "other_flags",
            "syn_type",
            "former_ranks",
            "problematic_synonyms",
            "synonyms",
        ]

    @property
    def valid_combination(self):
        return None if self.has_name is None else self.has_name.valid_combination

    def claim_problematic_synonym_statement(self, name, syn_type, error_str):
        if self.problematic_synonyms is None:
            self.problematic_synonyms = []
        blob = {"name": name, "syn_type": syn_type, "problem": error_str}
        self.problematic_synonyms.append(blob)

    def claim_type_material(self, type_str):
        n = self.name_attached_to_type_specimen
        try:
            n.claim_type_material(type_str)
        except:
            e = "could not find name that was type-material-based attached to TaxonConcept"
            self.claim_problematic_synonym_statement(type_str, "type material", e)

    @property
    def name_attached_to_type_specimen(self):
        return (
            None
            if self.has_name is None
            else self.has_name.name_attached_to_type_specimen
        )

    def _get_next_syn_id(self):
        ns = len(self.synonyms) if self.synonyms else 0
        try:
            trimmed_canon = ":".join(self.canonical_id.split(":")[2:])
        except:
            try:
                trimmed_canon = ":".join(self.canonical_id.split(":")[1:])
            except:
                trimmed_canon = self.canonical_id
        return "{}:syn{}".format(trimmed_canon, ns)

    @property
    def is_synonym(self):
        return self._is_synonym_of is not None

    def claim_is_synonym_of(self, valid):
        self._is_synonym_of = valid

    def claim_former_rank(self, rank):
        if self.former_ranks is None:
            self.former_ranks = []
        if rank not in self.former_ranks:
            self.former_ranks.append(rank)

    def claim_uninomial_synonym(self, name, syn_type, **kwargs):
        if "genus" in kwargs:
            return self._claim_syn_impl(name, syn_type, "genus", **kwargs)
        if "higher_group_name" not in kwargs:
            raise ValueError(
                "Expecting 'higher_group_name' or 'genus' kwarg "
                "in claim_uninomial_synonym"
            )
        return self._claim_syn_impl(name, syn_type, "clade", **kwargs)

    def claim_formerly_full_species(self, res, name, syn_type, **kwargs):
        self.claim_former_rank("species")
        return self.claim_binom_synonym(res, name, syn_type, **kwargs)

    def claim_formerly_subspecies(self, res, name, syn_type, **kwargs):
        self.claim_former_rank("subspecies")
        return self.claim_trinomial_synonym(res, name, syn_type, **kwargs)

    def claim_binom_synonym(self, res, name, syn_type, **kwargs):
        for expected in ("genus", "sp_epithet"):
            if expected not in kwargs:
                raise ValueError(
                    "Expecting '{}' kwarg in claim_binom_synonym".format(expected)
                )
        return self._claim_syn_impl(name, syn_type, "species", **kwargs)

    def claim_trinomial_synonym(self, res, name, syn_type, **kwargs):
        for expected in ("genus", "sp_epithet", "infra_epithet"):
            if expected not in kwargs:
                raise ValueError(
                    "Expecting '{}' kwarg in claim_trinomial_synonym".format(expected)
                )
        return self._claim_syn_impl(name, syn_type, "infraspecies", **kwargs)

    def _add_to_syn_list(self, syntc):
        if self.synonyms is None:
            self.synonyms = []
        self.synonyms.append(syntc)

    def _claim_syn_impl(self, name, syn_type, rank, **kwargs):
        tc = self.graph.add_taxon_concept(self._get_next_syn_id())
        self._add_to_syn_list(tc)
        tc.claim_rank(rank)
        tc.claim_is_synonym_of(self)
        if syn_type:
            tc.syn_type = syn_type
        from ..cmds.semanticize import semanticize_names

        semanticize_names(tc, name, kwargs, None)
        return tc
