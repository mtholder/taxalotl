#!/usr/bin/env python
from .graph_node import SemGraphNode


class NameSemNode(SemGraphNode):
    name_sem_nd_pred = ('name', 'authority')

    def __init__(self, sem_graph, id_minting_d, name):
        super(NameSemNode, self).__init__(sem_graph, id_minting_d)
        self._name = name
        self._authority = None

    @property
    def name(self):
        return self._name

    @property
    def predicates(self):
        return NameSemNode.name_sem_nd_pred

    def claim_authority(self, auth):
        if self._authority is None:
            self._authority = auth
        else:
            if not isinstance(self._authority, list):
                self._authority = [self._authority]
            self._authority.append(auth)

    @property
    def authority(self):
        if self._authority is None:
            return None
        if isinstance(self._authority, list):
            return [i.canonical_id for i in self._authority]
        return self._authority.canonical_id



class CombinationSemNode(NameSemNode):
    def __init__(self, sem_graph, id_minting_d, name):
        super(CombinationSemNode, self).__init__(sem_graph, id_minting_d, name)