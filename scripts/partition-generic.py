#!/usr/bin/env python
import sys
import os
def read_inp_taxonomy(fp, root_id):
    by_par_id, headers = {}, None
    root_line = None
    with open(fp, 'r', encoding='utf-8') as inp:
        par_id_ind, uid_ind = None, None
        for line in inp:
            fields = [i.strip() for i in line[:-1].split('\t|\t')]
            if par_id_ind is None:
                headers = fields
                par_id_ind = headers.index('parent_uid')
                uid_ind = headers.index('uid')
            else:
                par_id = fields[par_id_ind]
                by_par_id.setdefault(par_id, []).append(fields)
                if fields[uid_ind] == root_id:
                    assert root_line is None
                    root_line = fields
    all_ids = set()
    retained = {}
    if not root_line:
        return headers, retained, all_ids, None
    root_par_id = root_line[par_id_ind]
    retained[root_par_id] = [root_line]
    new_ids = [root_id]
    while new_ids:
        all_ids.update(new_ids)
        curr = list(new_ids)
        new_ids = []
        for cid in curr:
            children = by_par_id.get(cid)
            if children:
                retained[cid] = children
                new_ids.extend([i[uid_ind] for i in children])
    return headers, retained, all_ids, root_par_id

def read_inp_synonyms(fp, id_set):
    by_valid_id, headers = {}, None
    with open(fp, 'r', encoding='utf-8') as inp:
        uid_ind = None
        for line in inp:
            fields = [i.strip() for i in line[:-1].split('\t|\t')]
            if uid_ind is None:
                headers = fields
                uid_ind = headers.index('uid')
            else:
                uid = fields[uid_ind]
                if uid in id_set:
                    by_valid_id.setdefault(uid, []).append(fields)
    return headers, by_valid_id

def _write_line(out, fields):
    out.write('{}\n'.format('\t|\t'.join(fields)))

def write_taxonomy(fp, tax_headers, taxa, root_par_id):
    uid_ind = tax_headers.index('uid')
    with open(fp, 'w', encoding='utf-8') as outp:
        _write_line(outp, tax_headers)
        next_ids = [root_par_id]
        while next_ids:
            curr_ids = next_ids
            next_ids = []
            for i in curr_ids:
                child_list = taxa.get(i)
                if child_list:
                    for c in child_list:
                        _write_line(outp, c)
                        try:
                            uid = c[uid_ind]
                        except:
                            pass
                        else:
                            next_ids.append(uid)

def write_synonyms(fp, syn_headers, syns):
    with open(fp, 'w', encoding='utf-8') as outp:
        _write_line(outp, syn_headers)
        sk = list(syns.keys())
        sk.sort()
        for k in sk:
            for syn in syns[k]:
                _write_line(outp, syn)

def main():
    try:
        inp_dir, tag, root_id, out_dir = sys.argv[1:]
    except:
        sys.exit('Expected 4 arguments: input_dir, resource_tag, root_id, output_dir')
    tax_fn, syn_fn = 'taxonomy.tsv', 'synonyms.tsv'
    tax_blob = read_inp_taxonomy(os.path.join(inp_dir, tax_fn), root_id)
    tax_headers, tax_retained, tax_ids, root_par_id = tax_blob
    syn_blob = read_inp_synonyms(os.path.join(inp_dir, syn_fn), tax_ids)
    syn_headers, syn_by_valid_id = syn_blob
    write_taxonomy(os.path.join(out_dir, tax_fn), tax_headers, tax_retained, root_par_id)
    write_synonyms(os.path.join(out_dir, syn_fn), syn_headers, syn_by_valid_id)

if __name__ == '__main__':
    main()