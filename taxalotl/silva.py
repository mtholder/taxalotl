#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Previous silva handling code that served as a basis for this code was written by JAR and
#   Jessica Grant as a part of the reference_taxonomy and OToL efforts.
from __future__ import print_function

from peyotl import (assure_dir_exists,
                    get_logger)
import codecs
import os
from taxalotl.interim_taxonomy_struct import InterimTaxonomyData
from taxalotl.commands import unpack_resources
_LOG = get_logger(__name__)

def parse_silva_ids(fn):
    preferred = set()
    with codecs.open(fn, 'r', encoding='utf-8') as inp:
        for line in inp:
            ls = line.strip()
            if ls:
                preferred.add(ls)
    return preferred

def normalize_silva_taxonomy(source, destination, res_wrapper):
    assure_dir_exists(destination)
    itd = InterimTaxonomyData()
    id_list_id = res_wrapper.id_list
    taxalotl_config = res_wrapper.config
    id_list_res = taxalotl_config.get_terminalized_res_by_id(id_list_id, 'normalize silva')
    if not id_list_res.has_been_unpacked():
        unpack_resources(taxalotl_config, [id_list_id])
    expect_id_fp = os.path.join(id_list_res.unpacked_filepath, id_list_res.local_filename)
    if not os.path.isfile(expect_id_fp):
        raise ValueError("Silva ID file not found at: {}". format(expect_id_fp))
    expect_tax_fp = os.path.join(res_wrapper.unpacked_filepath, res_wrapper.local_filename)
    if not os.path.isfile(expect_tax_fp):
        raise ValueError("Silva taxon file not found at: {}". format(expect_tax_fp))

    preferred = parse_silva_ids(expect_id_fp)
    itd = InterimTaxonomyData()
    parse_silva_taxon_file(expect_tax_fp, preferred, itd)
    _LOG.info('{} taxonomy IDs read'.format(len(itd.to_par)))
    itd.write_to_dir(destination)

def gen_all_namepaths(path, name, prim_acc):
    an = []
    while path.endswith(';'):
        path = path[:-1]
    an.append(((path, name), prim_acc))
    ps = path.split(';')
    prev = ''
    for n, el in enumerate(ps):
        if not el:
            continue
        one_based = 1 + n
        np = (prev, el)
        prop_id = '{}/#{}'.format(prim_acc, one_based)
        an.append((np, prop_id))
        if prev:
            prev = '{};{}'.format(prev, el)
        else:
            prev = el
    return an

def parse_silva_taxon_file(expect_tax_fp, preferred_ids, itd):
    fung_pref = 'Eukaryota;Opisthokonta;Nucletmycea;Fungi;'
    animal_pref = 'Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);'
    plant_pref = 'Eukaryota;Archaeplastida;Chloroplastida;Charophyta;Phragmoplastophyta;Streptophyta;'
    mito_pref = 'Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria;'
    chloro_pref = 'Bacteria;Cyanobacteria;Chloroplast;'
    EUK = 'Eukaryota;'
    namepath_to_id_pair = {}
    with codecs.open(expect_tax_fp, 'r', encoding='utf-8') as inp:
        eh = 'primaryAccession\tstart\tstop\tpath\torganism_name\ttaxid\n'
        iinp = iter(inp)
        h = iinp.next()
        if h != eh:
            raise ValueError("Unexpected header: {}".format(h))
        for n, line in enumerate(iinp):
            ls = line.strip()
            if not ls:
                continue
            prim_acc, start, stop, path, name, tax_id = ls.split('\t')

            if n % 10000 == 0:
                _LOG.info("read taxon {} '{}' ...".format(n, name))
            if path.startswith(EUK):
                if (path.startswith(animal_pref)
                    or path.startswith(plant_pref)
                     or path.startswith(fung_pref)):
                    continue
            elif path.startswith('Bacteria;'):
                if path.startswith(mito_pref) or path.startswith(chloro_pref):
                    continue
            all_names = gen_all_namepaths(path, name, prim_acc)
            for np in all_names:
                #_LOG.info('np = {}'.format(np))
                assert not np[0][0].endswith(';')
                name_path, proposed_id = np
                stored = namepath_to_id_pair.setdefault(name_path, [None, None])
                if stored[0] is None:
                    if proposed_id in preferred_ids:
                        stored[0], stored[1] = proposed_id, proposed_id
                    elif stored[1] is None or proposed_id < stored[1]:
                        stored[1] = proposed_id
    for pid in namepath_to_id_pair.values():
        if pid[0] is None:
            pid[0] = pid[1]
    to_par = itd.to_par
    to_children = itd.to_children
    to_name = itd.to_name
    npk = namepath_to_id_pair.keys()
    npk.sort()
    to_par["0"] = None
    to_children["0"] = []
    to_name["0"] = "Life"
    itd.root_nodes.add("0")
    for name_path in npk:
        silva_id = namepath_to_id_pair[name_path][0]
        par_name = name_path[0]
        if par_name:
            if ';' in par_name:
                pnpl = par_name.split(';')
                pn = pnpl[-1]
                anc_part = ';'.join(pnpl[:-1])
                pnp = (anc_part, pn)
            else:
                pnp = ('', par_name)
            par_silva_id = namepath_to_id_pair[pnp][0]
        else:
            par_silva_id = "0"
        if silva_id in to_par:
            m = '{} remains mapped to ({}, {}) rather than ({}, {})'
            _LOG.warn(m.format(silva_id, to_par[silva_id], to_name[silva_id], par_silva_id, name_path[1]))
        else:
            to_par[silva_id] = par_silva_id
            to_children.setdefault(par_silva_id, []).append(silva_id)
            to_name[silva_id] = name_path[1]
    _LOG.info('{} SILVA ids stored'.format(len(to_name)))