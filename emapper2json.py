import sys
import json
from collections import defaultdict, Counter
from pymongo import MongoClient
from multiprocessing import Pool
import re

# #query_name     seed_eggNOG_ortholog    seed_ortholog_evalue    seed_ortholog_score     eggNOG OGs      narr_og_name    narr_og_cat     narr_og_desc    best_og_name    best_og_cat       best_og_desc    Preferred_name  GOs     EC      KEGG_ko KEGG_Pathway    KEGG_Module     KEGG_Reaction   KEGG_rclass     BRITE   KEGG_TC CAZy    BiGG_Reaction     PFAMs

col2name = {
            0: 'q', #query_name
            1: 'o',  #seed_eggNOG_ortholog
            2: 'ev', #seed_ortholog_evalue
            3: 'sc', #seed_ortholog_score
            4: 'ogs', #p_ogs VARCHAR,
            5: '?',#narr_og_name
            6: '?',#narr_og_cat
            7:'?',#narr_og_desc
            8:'best_og_name',#best_og_name
            9:'boc', #COG_category
            10:'?',  #Description
            11: 'pname', #p_name VARCHAR,
            12: 'gos', #p_go VARCHAR,
            13: 'ecs', #p_ec VARCHAR,
            14: 'kos', #p_ko VARCHAR,
            15: 'kpath', #p_kpath VARCHAR,
            16: 'kmods', #p_kmod VARCHAR,
            17: 'kreac', #kreact VARCHAR,
            18: 'krcls', #p_kclass VARCHAR,
            19: 'brite', #p_brite VARCHAR,
            20: 'ktc', #p_ktc VARCHAR,
            21: 'cazy', #p_cazy VARCHAR,
            22: 'bigg', #p_biggreact VARCHAR,
            23: 'pfam' #tax_scope VARCHAR,
        }

# get genes to report
hits_texas = set()
for line in open("../neighs/neigh_genes.txt"):
    hits_texas.add(line.rstrip())


# get annotations from mgv1 db
client = MongoClient('10.0.3.1')
db = client['mgv1']
col_emapper = db.emapper2

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

for i in chunks(list(hits_texas),100):
    matches = col_emapper.find({'q_g': {'$in': list(i)}})
    for m in matches:
        del m['_id']
        print (json.dumps(m))

# get annotations from current genomes
fname = sys.argv[1]
for line in open(fname):
    if line.startswith('#'):
        continue

    fields = list(map(str.strip, line.split('\t')))

    doc = {}
    for fn, f in enumerate(fields):
        if not f or f == '-':
            continue

        kname = col2name[fn]
        if kname in set(['q']):
            gene = f.split('@')[2]
            doc['q_g'] = gene

        elif kname in set(['ev', 'sc']):
            f = float(f)
        elif kname in set(['son', 'bon']):
            ogname, lvname = f.split('|')
            ogname, lvid = ogname.split('@')
            f = {'n': ogname, 'lvid': lvid, 'lvn': lvname}

        elif kname in set(['ogs']):
            og_array = []
            f = f.replace('dsDNA viruses, no RNA stage', 'dsDNAviruses-noRNAstage')
            for og in f.split(','):
                try:
                    ogname, level = og.split('|')
                except ValueError:
                    print(f, file=sys.stderr)
                    raise
                og_array.append(ogname)
            f = og_array

        elif kname in set(['kos']):
            f = [e.replace('ko:', '') for e in f.split(',')]

        elif kname in set(['sog', 'boc', 'gos', 'ecs', 'kpath', 'kmods', 'kreac', 'krcls', 'brite', 'ktc', 'cazy', 'bigg', 'pfam']):

            f = f.split(',')

        doc[kname] = f
    print(json.dumps(doc))
