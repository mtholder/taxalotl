#!/usr/bin/env python3

from sqlalchemy import and_
from sqlalchemy.orm import sessionmaker
from taxalotl.col_ingest.col_db_model import (AccordingTo,
                                              Authority,
                                              Base,
                                              CoLLink,
                                              Dataset, Description,
                                              Edge,
                                              Modified,
                                              NameFragment,
                                              Rank,
                                              SciNameUID, Synonym, SynonymStatus,
                                              Taxon, TaxonName,
                                              TaxonConcept, TaxonSuppInfo, TaxUUID,
                                              ValidNameStatus,
                                              )

from taxalotl.col_ingest.ingest_col_into_db import  get_engine

def taxon_from_name(session, name, allow_trailing_words=False):
    nw = [i.strip() for i in name.split()]
    if len(nw) < 1 or (len(nw) == 1 and not nw[0]):
        raise ValueError('expecting at least on word for "name", but got "{}"'.format(name))
    frag_lists = [session.query(NameFragment).filter(NameFragment.word == i).all() for i in nw]
    for n, resp in enumerate(frag_lists):
        if len(resp) != 1:
            raise ValueError('No name match for "{}"'.format(nw[n]))
    tnqa = [TaxonName.word_1_fk == frag_lists[0][0].pk]
    if allow_trailing_words:
        if len(frag_lists) > 1:
            tnqa.append(TaxonName.word_2_fk == frag_lists[1][0].pk)
            if len(frag_lists) > 2:
                tnqa.append(TaxonName.word_3_fk == frag_lists[2][0].pk)
                if len(frag_lists) > 3:
                    tnqa.append(TaxonName.word_4_fk == frag_lists[3][0].pk)
    else:
        tnqa.append(TaxonName.word_2_fk == frag_lists[1][0].pk if len(frag_lists) > 1 else None)
        tnqa.append(TaxonName.word_3_fk == frag_lists[2][0].pk if len(frag_lists) > 2 else None)
        tnqa.append(TaxonName.word_4_fk == frag_lists[3][0].pk if len(frag_lists) > 3 else None)
    return [session.query(NameFragment).get(i.word_1_fk).word for i in session.query(TaxonName).filter(and_(*tnqa)).all()]

if __name__ == '__main__':
    import sys
    from timeit import default_timer as timer
    engine = get_engine()
    session = sessionmaker(bind=engine)()
    start = timer()
    tlist = []
    for name in sys.argv[1:]:
        for taxon in taxon_from_name(session, name):
            tlist.append(taxon)
    end = timer()
    print('\n'.join([str(i) for i in tlist]))
    print(end - start, 'was the query time')
