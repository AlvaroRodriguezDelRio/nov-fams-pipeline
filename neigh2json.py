import json
import sys


target_files = sys.argv[1:]


# get json
for fname in target_files:
    for ln, line in enumerate(open(fname)):
        doc = {}
        contig, raw_genes = map(str.strip, line.split('\t'))
        raw_gene_set = set([i.split('|')[0] for i in raw_genes.split(',')])
        gene_array = []
        for index, gene_string in enumerate(raw_genes.split(',')):
            gene, start, end, strand = gene_string.split('|')
            if strand == '1':
                strand = '+'
            elif strand == '-1':
                strand = '-'
            d = {
                'g': gene.strip(),
                's': int(start),
                'e': int(end),
                'o': strand,
                'p': index,
            }
            gene_array.append(d)
            doc = {'c': contig,"genes":gene_array}
        print(json.dumps(doc))
