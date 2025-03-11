import sys
import re
from optparse import OptionParser

#args = sys.argv

parser = OptionParser()
parser.add_option("-e", "--extension", dest="ext", type="string",
                  help="fasta extension ('.faa','.fasta',etc)")
parser.add_option("-r", "--relpos", dest="pos", type="int",
                  help="relative position of genome name in file name")
parser.add_option("-t", "--taxonomy", dest="tax", type="string",
                  help="gtdb taxonomy file per genome")
parser.add_option("-p", "--path_fasta", dest="paths", type="string",
                  help="file with paths to fasta file per genome")

(options, args) = parser.parse_args()


# laod genome taxonomy
g2t = {}
for line in open(options.tax):
    f = list(map(str.strip,line.split('\t')))
    genome = f[0]
    tax = f[1]
    g2t[genome] = tax


for f in open(options.paths):
    f = f.rstrip()
    genome = f.split('/')[options.pos].replace(options.ext,'')
    tax = g2t[genome]
    d = tax.split(';')[0]
    p = tax.split(';')[1]
    if p == "p__":
        p = "NULL"
    for line in open(f):
        line = line.strip()
        if re.search(">",line):
            gene = line.replace('>','').split(' ')[0].split('\t')[0]
            gene_name = '@'.join(["HotSprings",genome,gene,d + '|' + p])
            print ('>' + gene_name)
        else:
            print(line.rstrip())
