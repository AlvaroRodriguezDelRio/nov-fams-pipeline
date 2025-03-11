import sys
import json
from collections import defaultdict, Counter
import re
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-e", "--emapper_file", dest="emapper_file", type="string",
                  help="functional annotation file (eggnog-mapper output)")
(options, args) = parser.parse_args()

col2name = {
            0: 'q', #query_name
            1: 'o',  #seed_eggNOG_ortholog
            2: 'ev', #seed_ortholog_evalue
            3: 'sc', #seed_ortholog_score
            4: 'ogs', #p_ogs VARCHAR,
            5: '?',#max_annot_lvl
            6:'boc', #COG_category
            7:'?',  #Description
            8: 'pname', #p_name VARCHAR,
            9: 'gos', #p_go VARCHAR,
            10: 'ecs', #p_ec VARCHAR,
            11: 'kos', #p_ko VARCHAR,
            12: 'kpath', #p_kpath VARCHAR,
            13: 'kmods', #p_kmod VARCHAR,
            14: 'kreac', #kreact VARCHAR,
            15: 'krcls', #p_kclass VARCHAR,
            16: 'brite', #p_brite VARCHAR,
            17: 'ktc', #p_ktc VARCHAR,
            18: 'cazy', #p_cazy VARCHAR,
            19: 'bigg', #p_biggreact VARCHAR,
            20: 'pfam' #tax_scope VARCHAR,
        }


# get annotations from current genomes
fname = options.emapper_file
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
