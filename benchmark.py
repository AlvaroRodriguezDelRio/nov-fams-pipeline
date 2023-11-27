import sys
from collections import defaultdict,Counter
from itertools import groupby
import itertools
import json

kpaths_val = (set(i.rstrip().split('\t')[0] for i in open('KEGG_pathways.tab')))
fams_more_2 = set([i.rstrip() for i in open('microbial_genomes-v1.clustering.folded.parsed.no_dups.fams_more_2_members.txt')])

# get gene families whose all members are annotated to the same term by homology
fam2annots = defaultdict(lambda:defaultdict(lambda:set()))
for line in open("KEGG pathway annotation per gene family.tab"):
    fam,db,ngenes,annot,n_annot,desc = list(map(str.strip,line.split('\t')))
    if int (ngenes) == int(n_annot) and int(ngenes) > 2 and fam in fams_more_2:
        fam2annots[fam][db].add(annot)

# combination of genomic context parameters to explore
v_cons  = [0.9,0.8,0.5,0.3]
h_cons = [1,2,3,4]
pos = [[1],[1,2],[1,2,3]]
op_strand = [0,0.2,1]
dist = [0,0.2,1]

# compare KEGG patways assigned by homology to predictions by genomic context
with open("score_per_pos.json") as jfile:
    for line in jfile:
        line = json.loads(line)
        fam = line['fam']

        if fam not in fams_more_2 or fam not in fam2annots:
            continue
          
        # iterate through KEGG pathway kpath vertical conservation values, per position
        for db,annots in line.items():
            if db in ['kpath']:
                preds = defaultdict(lambda:Counter())
                for annot in annots:
                    if annot['n'] in kpaths_val:

                        # iterate through conservation combinations that we want to contemplate
                        for comb in itertools.product(v_cons,pos,op_strand,dist):
                            if float(annot['score']) > comb[0] and float(annot['mean_num_in_opposite_strand']) <= comb[2] and float(annot['mean_num_pos_opposite_strand_between']) <= comb[2] and float(annot['num_h_dis']) <= comb[3] and abs(annot['pos']) in comb[1]:
                                comb_name = '_'.join(map(str,[comb[0],'-'.join(map(str,comb[1])),comb[2],comb[3]]))
                                preds[comb_name][annot['n']] += 1

                # compare homology and neigh predictions for this family, for different horizontal conservation scores
                for i in h_cons:
                    for c in preds:
                        f_preds = set()
                        for annot in preds[c]:
                            if preds[c][annot] >= i:
                                f_preds.add(annot)

                        if f_preds:
                            coincidences = f_preds.intersection(fam2annots[fam][db])
                            incorrect = f_preds - fam2annots[fam][db]
                            missing = fam2annots[fam][db] - f_preds
                            print ('\t'.join([fam,db,c,str(i),str(len(list(coincidences))),','.join(list(coincidences)),str(len(list(incorrect))),','.join(list(incorrect)),str(len(list(missing))),','.join(list(missing))]))
