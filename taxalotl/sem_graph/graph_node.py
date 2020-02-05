#!/usr/bin/env python

class SemGraphNode(object):
    def __init__(self, sem_graph, id_minting_d):
        self.canonical_id = sem_graph.register_obj(id_minting_d, self)
        self.graph = sem_graph

    def as_dict(self):
        d = {}
        for att in self.predicates:
            val = getattr(self, att, [])
            if val:
                if isinstance(val, list) or isinstance(val, tuple):
                    d[att] = [_serialize_triple_object(i) for i in val]
                else:
                    d[att] = _serialize_triple_object(val)
        return d

    @property
    def predicates(self):
        return []


def _serialize_triple_object(o):
    return o.canonical_id if isinstance(o, SemGraphNode) else o


class AuthoritySemNode(SemGraphNode):
    auth_sem_nd_pred = ('authors', 'year', 'taxon_concepts')

    def __init__(self, sem_graph, id_minting_d, authors, year, tax_con_sem_node):
        super(AuthoritySemNode, self).__init__(sem_graph, id_minting_d)
        self._authors = authors
        self._year = year
        self.taxon_concept_set = set([tax_con_sem_node])

    @property
    def taxon_concepts(self):
        return [i.canonical_id for i in self.taxon_concept_set]

    @property
    def year(self):
        return self._year

    @property
    def authors(self):
        return self._authors

    @property
    def predicates(self):
        return AuthoritySemNode.auth_sem_nd_pred

