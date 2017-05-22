#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
# Previous silva handling code that served as a basis for this code was written by JAR and
#   Jessica Grant as a part of the reference_taxonomy and OToL efforts.
from __future__ import print_function

from peyotl import (assure_dir_exists,
                    get_logger,
                    read_as_json,
                    write_as_json)
import codecs
import os
from taxalotl.interim_taxonomy_struct import InterimTaxonomyData
from taxalotl.partitions import do_partition, GEN_MAPPING_FILENAME
from taxalotl.commands import unpack_resources
from taxalotl.ott import _partition_ott_by_root_id

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
    depends_on = res_wrapper.depends_on
    taxalotl_config = res_wrapper.config
    expect_id_fp, ncbi_mapping_fp = None, None
    for dep_id in depends_on:
        dep_res = taxalotl_config.get_terminalized_res_by_id(dep_id, 'normalize silva')
        if not dep_res.has_been_unpacked():
            unpack_resources(taxalotl_config, [dep_id])
        dep_fp = os.path.join(dep_res.unpacked_filepath, dep_res.local_filename)
        if not os.path.isfile(dep_fp):
            raise ValueError("Silva ID file not found at: {}". format(dep_fp))
        if dep_res.schema.lower() == 'id list':
            expect_id_fp = dep_fp
        elif dep_res.schema.lower() == 'silva taxmap':
            ncbi_mapping_fp = dep_fp
        else:
            raise ValueError('unrecognized dependency schema {}'.format(dep_res.schema))
    if expect_id_fp is None:
        raise ValueError('ID list dependency not found.')
    if ncbi_mapping_fp is None:
        raise ValueError('NCBI mapping dependency not found.')
    expect_tax_fp = os.path.join(res_wrapper.unpacked_filepath, res_wrapper.local_filename)
    if not os.path.isfile(expect_tax_fp):
        raise ValueError("Silva taxon file not found at: {}". format(expect_tax_fp))
    acc_to_trim = parse_acc_to_trim_from_ncbi(ncbi_mapping_fp)
    preferred = parse_silva_ids(expect_id_fp)
    itd = InterimTaxonomyData()
    part_name_to_silva_id = parse_silva_taxon_file(expect_tax_fp, preferred, acc_to_trim, itd)
    _LOG.info('{} taxonomy IDs read'.format(len(itd.to_par)))
    itd.write_to_dir(destination)
    mapping_file = os.path.join(destination, GEN_MAPPING_FILENAME)
    write_as_json(part_name_to_silva_id, mapping_file, indent=2, separators=(',', ': '))

def parse_acc_to_trim_from_ncbi(tax_fp):

    trimmed_pref = ('root;cellular organisms;Eukaryota;Opisthokonta;Fungi;',
                    'root;cellular organisms;Eukaryota;Opisthokonta;Metazoa;',
                    'root;cellular organisms;Eukaryota;Viridiplantae;',)
    EUK = 'root;cellular organisms;Eukaryota;'
    to_trim = set()
    with codecs.open(tax_fp, 'r', encoding='utf-8') as inp:
        for n, line in enumerate(inp):
            ls = line.strip()
            if not ls:
                continue
            prim_acc, start, stop, path, name = ls.split('\t')
            if n % 10000 == 0:
                _LOG.info("scanned taxon {} '{}' ...".format(n, name))
            if path.startswith(EUK):
                for pref in trimmed_pref:
                    if path.startswith(pref):
                        to_trim.add(prim_acc)
                        break
    return to_trim

def gen_all_namepaths(path, name, prim_acc):
    an = []
    while path.endswith(';'):
        path = path[:-1]
    if name:
        if '(' in name:
            name = name.split('(')[0].strip()
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

def parse_silva_taxon_file(expect_tax_fp, preferred_ids, acc_to_trim, itd):
    fung_pref = 'Eukaryota;Opisthokonta;Nucletmycea;Fungi;'
    animal_pref = 'Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);'
    plant_pref = 'Eukaryota;Archaeplastida;Chloroplastida;Charophyta;Phragmoplastophyta;Streptophyta;'
    mito_pref = 'Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria;'
    chloro_pref = 'Bacteria;Cyanobacteria;Chloroplast;'
    trim_pref = (fung_pref, animal_pref, plant_pref, mito_pref, chloro_pref)

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
            if prim_acc in acc_to_trim:
                tpath = None
                for p in trim_pref:
                    if path.startswith(p):
                        tpath = p
                        break
                if tpath is None:
                    _LOG.info('deleting to untrimmable {}'.format(prim_acc))
                    continue
                else:
                    all_names = gen_all_namepaths(path, '', prim_acc)
            else:
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
    part_map = {
        'Archaeplastida': ('Eukaryota', 'Archaeplastida'),
        'Chloroplastida': ('Eukaryota;Archaeplastida', 'Chloroplastida'),
        'Glaucophyta': ('Eukaryota;Archaeplastida', 'Glaucophyta'),
        'Rhodophyta': ('Eukaryota;Archaeplastida', 'Rhodophyceae'),
        'Archaeplastida': ('Eukaryota', 'Archaeplastida'),
        'Haptophyta': ('Eukaryota', 'Haptophyta'),
        'Eukaryota':  ('', 'Eukaryota'),
        'Archaea': ('', 'Archaea'),
        'Bacteria': ('', 'Bacteria'),
    }
    for pid in namepath_to_id_pair.values():
        if pid[0] is None:
            pid[0] = pid[1]
    part_map_to_namepath = {
        'Archaeplastida': ('Eukaryota', 'Archaeplastida'),
        'Chloroplastida': ('Eukaryota;Archaeplastida', 'Chloroplastida'),
        'Glaucophyta': ('Eukaryota;Archaeplastida', 'Glaucophyta'),
        'Rhodophyta': ('Eukaryota;Archaeplastida', 'Rhodophyceae'),
        'Archaeplastida': ('Eukaryota', 'Archaeplastida'),
        'Haptophyta': ('Eukaryota', 'Haptophyta'),
        'Eukaryota': ('', 'Eukaryota'),
        'Archaea': ('', 'Archaea'),
        'Bacteria': ('', 'Bacteria'),
        'SAR': ('Eukaryota', 'SAR'),
    }
    part_name_to_silva_id = {}
    for part_name, path_name in part_map_to_namepath.items():
        part_name_to_silva_id[part_name] = [namepath_to_id_pair[path_name][0]]
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
    return part_name_to_silva_id

def partition_silva(res_wrapper, part_name, part_keys, par_frag):
    silva_map = read_as_json(os.path.join(res_wrapper.normalized_filepath, GEN_MAPPING_FILENAME))
    do_partition(res_wrapper,
                 part_name,
                 part_keys,
                 par_frag,
                 master_map=silva_map,
                 parse_and_partition_fn=_partition_ott_by_root_id)


