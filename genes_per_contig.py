import sys
import re
from collections import defaultdict

gene2scaff = defaultdict(lambda:defaultdict(lambda:[]))
contig2genome = {}
for file_n in open(sys.argv[1]):
    file_n = file_n.rstrip()
    genome = file_n.split('\t')[-1].replace('.gff','')
    for line in open(file_n):
        if re.search("ID=",line):
                contig = line.split('\t')[0]
                contig2genome[contig] = genome
                gene = line.split('\t')[8].split(";")[0].replace("ID=","")
                beg = line.split('\t')[3]
                end = line.split("\t")[4]
                strand = line.split("\t")[6]
                gene2scaff[contig][int(beg)].append([gene,int(end),strand])


for contig in gene2scaff:
        scaff_starts_sorted = sorted(gene2scaff[contig].keys())
        neighs_array = []
        for gene_pos in scaff_starts_sorted:
                for elem in gene2scaff[contig][gene_pos]:
                        target_gene = elem
                        target_gene_name = target_gene[0] + "|" + str(gene_pos) + "|" + str(target_gene[1]) + "|" + str(target_gene[2])
                neighs_array.append(target_gene_name)
        print ('\t'.join([contig,','.join(neighs_array)]))
