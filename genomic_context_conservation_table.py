import json
import sys
from itertools import groupby
from collections import defaultdict

def readlines(file_):
    for line in open(file_):
        yield list(map(str.strip,line.split('\t')))


def calculate_number_neg_between(pos,neg_array,positions_annot):
    number_neg = 0
    pos = int(pos)
    positions_annot = [x for x in positions_annot.split(',')]

    if int(pos)>0:
        for i in range(1,pos):
            for j in neg_array[str(i)]:
                if str(j) in positions_annot:
                    number_neg += 1
    elif pos <0:
        for i in range(pos+1,0):
            for j in neg_array[str(i)]:
                if str(j) in positions_annot:
                    number_neg += 1
    return number_neg

with open(sys.argv[1]) as file1:
    for k, g in groupby(file1, key=lambda x:x.split('\t')[0]):
        neg_strand_pos = defaultdict(lambda:set())
        more_x_nucleotide_distance = defaultdict(lambda:set())
        lines = []
        for line in g:
            if len(line.split('\t')) < 3 or line == '':
                continue

            fam,pos,db,annot,number_genes_with_annot,number_genes_per_pos,nseqs_fam,positions_annot,positions_annot_contrary,espace,num_contrary_strand = list(map(str.strip,line.split('\t')))

            # gather contrary strand coordinates
            if int (num_contrary_strand) > 0:
                for i in positions_annot_contrary.split(','):
                    neg_strand_pos[pos].add(int(i))


            # load genes far away from the gene closer to the anchor
            for i,dist in enumerate(espace.split(',')):
                if int(dist) > 100:
                    n_gene_dist = positions_annot.split(',')[i]
                    more_x_nucleotide_distance[int(pos)].add(int(n_gene_dist))


            if db != 'unknown':
                lines.append(line)


        # accumulate high distances until the anchor
        # if the -3 is far from -2, and -2 is far from -1, then -3 is 'double distant'
        neg_pos = [-3,-2,-1]
        for x in neg_pos:
            for y in neg_pos:
                if x < y: # add all high distance positions before anchor and gene
                    more_x_nucleotide_distance[x] = list(more_x_nucleotide_distance[x]) # convert to array for including more distant distanced (len of the array used)
                    for z in more_x_nucleotide_distance[y]:
                        more_x_nucleotide_distance[x].append(z)

        pos_pos = [3,2,1]
        for x in pos_pos:
            for y in pos_pos:
                if x > y:  # add all high distance positions before anchor and gene
                    more_x_nucleotide_distance[x] = list(more_x_nucleotide_distance[x]) # convert to array for including more distant distanced (len of the array used)
                    for z in more_x_nucleotide_distance[y]:
                        more_x_nucleotide_distance[x].append(z)

        added_cogs = set() # avoid reporting for the same family and por the same cogs several times
        fam2info = defaultdict(lambda:[])
        for line in lines:
            fam,pos,db,annot,number_genes_with_annot,number_genes_per_pos,nseqs_fam,positions_annot,positions_annot_contrary,espace,num_contrary_strand = list(map(str.strip,line.split('\t')))
            fam2info['fam'] = fam
            if annot.split('@')[0]+pos not in added_cogs:
                added_cogs.add(annot.split('@')[0] + pos)
                score = int(number_genes_with_annot)/int(nseqs_fam)
                number_neg_between = calculate_number_neg_between(pos,neg_strand_pos,positions_annot)
                number_genes_with_annot = len(positions_annot.split(','))
                number_genes_high_distance = len([x for x in more_x_nucleotide_distance[int(pos)] if str(x) in positions_annot.split(',')])
                desc = ''
                if annot.split('@')[0] in kegg2annot:
                    desc = kegg2annot[annot.split('@')[0]]
                print ('\t'.join(list(map(str,[fam,db,annot.split('@')[0],pos,score,int(num_contrary_strand)/int(number_genes_with_annot),int(number_neg_between)/int(number_genes_with_annot),int(number_genes_high_distance)/int(number_genes_with_annot),desc]))))
                cog_dict = {'n':annot.split('@')[0],'score':float(score),"mean_num_in_opposite_strand":float(int(num_contrary_strand)/int(number_genes_with_annot)),"mean_num_pos_opposite_strand_between":float(int(number_neg_between)/int(number_genes_with_annot)),"pos":int(pos),"num_h_dis":float(int(number_genes_high_distance)/int(number_genes_with_annot))}
                fam2info[db].append(cog_dict)
