from ete3 import SeqGroup
import sys
from multiprocessing import Pool


def iter(paths):
        with open(paths) as table:
                for line in table:
                        line = line[:-1]
                        yield (line)


def get_domain_len(line):
        alg = SeqGroup(line)
        alg_len = 0
        name2sequence = dict()
        number_seqs = 0
        for sequence in alg:
                name = sequence[0]
                alg_len = len(alg.get_seq(name))
                name2sequence[name] = alg.get_seq(name)
                number_seqs += 1
        number_no_gaps = dict()
        for index in range(0,alg_len):
                n_no_gaps = 0
                for gene_name in name2sequence:
                        if name2sequence[gene_name][index] != '-':
                                n_no_gaps += 1
                number_no_gaps[index] = n_no_gaps

        consecutive = 0
        max_consecutive = 0
        for i in number_no_gaps:
                #print (consecutive, max_consecutive)
                if float(number_no_gaps[i]) > float(number_seqs)*0.8:
                        consecutive += 1
                else:
                        consecutive = 0
                if consecutive > max_consecutive:
                        max_consecutive = (consecutive)


        cluster_name = line.split('/')[-1]
        return (cluster_name,float(max_consecutive)/float(alg_len), max_consecutive)



algfile_paths = sys.argv[1]

pool=Pool(1)
out_file = open("alg_quality.number_conseq_res_domain.tab","w")
for result in pool.imap(get_domain_len,iter(algfile_paths)):
        out_file.write('\t'.join(map(str,result)))
        out_file.write('\n')
