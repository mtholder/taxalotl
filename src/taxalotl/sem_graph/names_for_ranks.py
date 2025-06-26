#!/usr/bin/env python
from .graph_node import SemGraphNode
from .name import NameSemNode

class TypeSpecimen(SemGraphNode):
    sem_nd_pred = ('specimen_codes', 'name_for_type', 'valid_name')

    @property
    def predicates(self):
        return TypeSpecimen.sem_nd_pred

    def __init__(self, sem_graph, id_minting_d, spec_code, name_for_type, valid_name):
        id_minting_d['class_tag'] = 'type'
        super(TypeSpecimen, self).__init__(sem_graph, id_minting_d)
        self.name_for_type = name_for_type
        self.specimen_codes = spec_code
        self.valid_name = valid_name

class GenusGroupSemNode(NameSemNode):
    def __init__(self, sem_graph, id_minting_d, name):
        super(GenusGroupSemNode, self).__init__(sem_graph, id_minting_d, name)
        self.contained = []


class SpeciesGroupSemNode(NameSemNode):
    sp_grp_name_sem_nd_pred = tuple(list(NameSemNode.name_sem_nd_pred) + ['authority'])

    def __init__(self, sem_graph, id_minting_d, name):
        super(SpeciesGroupSemNode, self).__init__(sem_graph, id_minting_d, name)
        self._authority = None
        self.contained = []


class HigherGroupSemNode(NameSemNode):
    def __init__(self, sem_graph, id_minting_d, name):
        super(HigherGroupSemNode, self).__init__(sem_graph, id_minting_d, name)


class SpecimenCodeSemNode(NameSemNode):
    def __init__(self, sem_graph, id_minting_d, name):
        super(SpecimenCodeSemNode, self).__init__(sem_graph, id_minting_d, name)
