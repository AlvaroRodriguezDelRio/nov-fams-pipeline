import sys
import re
from collections import defaultdict
from optparse import OptionParser
import json

parser = OptionParser()
parser.add_option("-e", "--extension", dest="ext", type="string",
                  help="fasta extension ('.faa','.fasta',etc)")
parser.add_option("-r", "--relpos", dest="pos", type="int",
                  help="relative position of genome name in file name")
parser.add_option("-p", "--path_fasta", dest="paths", type="string",
                  help="file with paths to fasta file per genome")

(options, args) = parser.parse_args()


gene2scaff = defaultdict(lambda:defaultdict(lambda:[]))
contig2genome = {}
for file_n in open(options.paths):
    file_n = file_n.rstrip()
    genome = file_n.split('/')[options.pos].replace(options.ext,'')
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
    doc = {}
    gene_array = []
    for i,gene_pos in enumerate(scaff_starts_sorted):
        for gene_info in gene2scaff[contig][gene_pos]:
            #target_gene_name = target_gene[0] + "|" + str(gene_pos) + "|" + str(target_gene[1]) + "|" + str(target_gene[2])
            gene = gene_info[0]
            start =  str(gene_pos)
            end = str(gene_info[1])
            strand = str(gene_info[2])

            if strand == '1':
                strand = '+'
            elif strand == '-1':
                strand = '-'
            d = {
                'g': gene.strip(),
                's': int(start),
                'e': int(end),
                'o': strand,
                'p': i + 1,
            }
            gene_array.append(d)
    doc = {'c': contig,"genes":gene_array}
    print(json.dumps(doc))
