#!/usr/bin/env python
from .name import NameSemNode


class VerbatimSemNode(NameSemNode):
    extra_pred = ('combination',
                  'genus_name',
                  'higher_group_name',
                  'infra_epithets',
                  'normalized',
                  'sp_epithet',
                  'specimen_codes',
                  'subgenus_names',
                  )

    name_sem_nd_pred = tuple(list(NameSemNode.name_sem_nd_pred) + list(extra_pred))

    def __init__(self, sem_graph, id_minting_d, name):
        super(VerbatimSemNode, self).__init__(sem_graph, id_minting_d, name)
        self.combination = None
        self.genus_name = None
        self.higher_group_name = None
        self.infra_epithets = None
        self.normalized = None
        self.sp_epithet = None
        self.specimen_codes = None
        self.subgenus_names = None

    @property
    def canonical_name(self):
        if self.combination is not None:
            return self.combination
        if self.subgenus_names is not None:
            return self.subgenus_names
        if self.genus_name is not None:
            return self.genus_name
        if self.higher_group_name is not None:
            return self.higher_group_name
        assert self.name
        return self

    @property
    def predicates(self):
        return VerbatimSemNode.name_sem_nd_pred

    def claim_higher_group_name(self, n):
        assert self.higher_group_name is None or self.higher_group_name.name == n.name
        if self.higher_group_name and self.higher_group_name.name == n.name:
            return
        self.higher_group_name = n

    def claim_combination(self, n):
        assert self.combination is None or self.combination is n
        self.combination = n

    def claim_normalized(self, n):
        assert self.normalized is None or self.normalized is n
        self.normalized = n

    def claim_genus(self, n):
        assert self.genus_name is None or self.genus_name == n
        self.genus_name = n

    def claim_subgenus(self, n):
        if self.subgenus_names is None:
            self.subgenus_names = [n]
        else:
            self.subgenus_names.append(n)

    def claim_sp_epithet(self, n):
        assert self.sp_epithet is None or self.sp_epithet.name == n.name
        if self.sp_epithet and self.sp_epithet.name == n.name:
            return
        self.sp_epithet = n

    def claim_infra_epithet(self, n):
        if self.infra_epithets is None:
            self.infra_epithets = [n]
        else:
            self.infra_epithets.append(n)

    def claim_specimen_code(self, n):
        if self.specimen_codes is None:
            self.specimen_codes = [n]
        else:
            self.specimen_codes.append(n)

    @property
    def valid_combination(self):
        return None if self.combination is None else self.combination.name

    @property
    def most_terminal_name(self):
        ie = self.most_terminal_infra_epithet
        if ie is not None:
            return ie
        for n in [self.sp_epithet, self.subgenus_names, self.genus_name, self.higher_group_name]:
            if n:
                return n
        return None

    @property
    def most_terminal_infra_epithet(self):
        if self.infra_epithets:
            if len(self.infra_epithets) > 1:
                x = [(len(i.name), i.name, i) for i in self.infra_epithets]
                x.sort()
                return x[-1][-1]
            return self.infra_epithets[0]
        return None

    @property
    def name_attached_to_type_specimen(self):
        ie = self.most_terminal_infra_epithet
        if ie is not None:
            return ie
        return self.sp_epithet

    def add_normalized(self, name_str):
        n = self.graph._add_normalized(self, name_str)
        self.claim_normalized(n)
        return n

    def add_combination(self, name_str):
        n = self.graph._add_combination(self, name_str)
        self.claim_combination(n)
        return n

    def add_genus(self, name_str):
        n = self.graph._add_genus(self, name_str)
        self.claim_genus(n)
        return n

    def add_subgenus(self, name_str):
        n = self.graph._add_subgenus(self, name_str)
        self.claim_subgenus(n)
        return n

    def add_sp_epithet(self, name_str, genus_sem_node, avoid_dup=True):
        n = self.graph._add_sp_epithet(self, name_str, genus_sem_node, avoid_dup=avoid_dup)
        self.claim_sp_epithet(n)
        return n

    def add_infra_epithet(self, name_str, epi_sem_node):
        n = self.graph._add_infra_epithet(self, name_str, epi_sem_node)
        self.claim_infra_epithet(n)
        return n

    def add_higher_group_name(self, name_str):
        n = self.graph._add_higher_group_name(self, name_str)
        self.claim_higher_group_name(n)
        return n

    def add_specimen_code(self, name_str):
        n = self.graph._add_specimen_code(self, name_str)
        self.claim_specimen_code(n)
        return n
