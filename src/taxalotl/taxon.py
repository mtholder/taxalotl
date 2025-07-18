#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from .taxonomic_ranks import _RANK_TO_SORTING_NUMBER


# end taxonomic ranks


class Taxon(object):
    _DATT = ("id", "par_id", "name", "rank", "src_dict", "flags", "uniqname")

    def __init__(self, line=None, line_num="<unknown>", line_parser=None, d=None):
        self.id, self.par_id, self.name, self.rank = None, None, None, None
        self.src_dict, self.flags, self.uniqname = None, None, None
        self.children_refs = None
        self._synonyms = None
        if d is not None:
            self.from_serializable_dict(d)
        else:
            self.line_num = line_num
            self.line = line
            if line_parser is None:
                from .ott_schema import full_ott_line_parser

                line_parser = full_ott_line_parser
            line_parser(self, line)

    def formatted_src_dict(self):
        if not self.src_dict:
            return
        sortable_list = [(k, list(v)) for k, v in self.src_dict.items()]
        sortable_list.sort()
        strs = []
        for el in sortable_list:
            el[1].sort()
            strs.extend(["{}:{}".format(el[0], i) for i in el[1]])
        return ",".join(strs)

    def child_id_dict(self):
        return {c.id: c for c in self.children_refs}

    def get_flag_str(self):
        sf = self.sorted_flags
        return ",".join(sf) if sf else ""

    def terse_descrip(self):
        m = '"{}" [{}]{}\n'
        suff = ""
        r = self.rank
        if r:
            suff += " {}".format(r)
        f = " flags={}".format(self.get_flag_str())
        if f:
            suff += f
        return m.format(self.name, self.id, suff)

    @property
    def sorted_flags(self):
        if not self.flags:
            return []
        tmp = list(self.flags)
        tmp.sort()
        return tmp

    @property
    def synonyms(self):
        if hasattr(self, "_synonyms"):
            return self._synonyms
        return set()

    @synonyms.setter
    def synonyms(self, x):
        self._synonyms = x

    def add_synonym_id(self, syn_id):
        if self._synonyms is None:
            self._synonyms = set()
        self._synonyms.add(syn_id)

    def flag_as_hybrid(self):
        self._add_flag("hybrid")

    def flag_as_incertae_sedis(self):
        self._add_flag("incertae_sedis")

    def _add_flag(self, f):
        if not self.flags:
            self.flags = {f}
        else:
            self.flags.add(f)

    def rank_sorting_number(self):
        if (self.rank is None) or self.rank.startswith("no rank"):
            return None
        return _RANK_TO_SORTING_NUMBER[self.rank]

    def __str__(self):
        s = "rank={}".format(self.rank) if self.rank else ""
        m = 'Taxon: "{}" (id={} | par={} | {})'
        return m.format(self.name, self.id, self.par_id, s)

    def __repr__(self):
        return "Taxon(d={})".format(self.to_serializable_dict())

    @property
    def name_that_is_unique(self):
        return self.uniqname if self.uniqname else self.name

    def from_serializable_dict(self, d):
        for k, v in d.items():
            if k == "flags":
                v = set(v)
            elif k == "src_dict":
                v = {sk: set(sv) for sk, sv in v.items()}
            setattr(self, k, v)

    def to_serializable_dict(self):
        d = {}
        for k in Taxon._DATT:
            v = getattr(self, k, None)
            if v is not None:
                if k == "id" or k == "par_id":
                    d[k] = v
                elif v:
                    if k == "flags":
                        if len(v):
                            d[k] = list(v)
                            d[k].sort()
                    elif k == "src_dict":
                        ds = {}
                        for sk, ss in v.items():
                            sl = list(ss)
                            sl.sort()
                            ds[sk] = sl
                        d[k] = ds
                    else:
                        d[k] = v
        return d

    def get(self, key, default=None):
        return getattr(self, key, default)

    def __getitem__(self, item):
        return self.__dict__[item]
