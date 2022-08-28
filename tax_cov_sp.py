import sys
from collections import defaultdict,Counter

lin2num = Counter()
genome2t = {}
for f in ['/data/jhc/cold/MAGs/combined/genome2taxonomy.r202.tab','../data/genome2taxonomy.custom.tab']:
    for line in open(f):
        genome,tax = list(map(str.strip,line.split('\t')))
        genome2t[genome] = tax
        genome2t[genome.split('.')[0]] = tax
        t_comb = []
        # number of genomes per lin (for calculating coverage)
        for t in tax.split(';'):
            t_comb.append(t)
            lin2num[';'.join(t_comb)] += 1

fam2genomes = defaultdict(lambda:set())
for line in open(sys.argv[1]):
    fam,n,mems = list(map(str.strip,line.split('\t')))
    n = len(mems.split(','))
    for m in mems.split(','):
        genome = m.split('@')[1].replace('.oneLine','')
        fam2genomes[fam].add(genome)

    fam_counter = Counter()
    for genome in fam2genomes[fam]:
        if genome not in genome2t:
            genome = 'GB_' + genome
        tax = genome2t[genome]
        t_comb = []
        # number of genomes per lin (for calculating coverage)
        for t in tax.split(';'):
            t_comb.append(t)
            fam_counter[';'.join(t_comb)] += 1

    for lin, number in fam_counter.items():
        sp = str(number / len(list(fam2genomes[fam])))
        cov = str(number / lin2num[lin])
        print ('\t'.join([fam,lin,sp,cov,str(lin2num[lin])]))
