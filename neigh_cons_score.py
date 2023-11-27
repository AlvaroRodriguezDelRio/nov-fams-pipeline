import sys
import json
from collections import defaultdict, Counter
from pymongo import MongoClient
from multiprocessing import Pool
import re


def process_fam(line):
    client = MongoClient('CLIENT CODE')
    db = client['DB_NAME']
    col_emapper = db.emapper2
    col_neighs = db.neighs
    col_cards = db.CARD

    def get_distances(minicontig,pos_anchor,mini_contig_genes,anchor_strand):
        distances = defaultdict(lambda:{})
        prev_end = 0
        prev_start = 0
        for i,gene in enumerate(mini_contig_genes):
            rel_pos = i - pos_anchor
            if anchor_strand[0] == '+':
                start = min([n['s'] for n in mini_contig if n['g'] == gene][0],[n['e'] for n in mini_contig if n['g'] == gene][0])
                end = max([n['s'] for n in mini_contig if n['g'] == gene][0],[n['e'] for n in mini_contig if n['g'] == gene][0])
                if prev_end != 0:
                    dis = start - prev_end
                    distances[rel_pos-1][rel_pos] = dis
                prev_end = end
                prev_start = start
            else:
                start = max([n['s'] for n in mini_contig if n['g'] == gene][0],[n['e'] for n in mini_contig if n['g'] == gene][0])
                end = min([n['s'] for n in mini_contig if n['g'] == gene][0],[n['e'] for n in mini_contig if n['g'] == gene][0])
                if prev_end != 0:
                    dis =  prev_end - start
                    distances[rel_pos-1][rel_pos] = dis
                prev_end = end

        return distances

    def get_emapper_annotations(names):
        matches = col_emapper.find({'q_g': {'$in': names} })
        gene2annot = defaultdict(dict)
        for m in matches:
            ogs_by_level = {}
            del m['_id']
            gene2annot[m['q_g']] = m
        return gene2annot


    def get_cards(names):
        matches = col_cards.find({'q_g': {'$in': names} })
        gene2card = defaultdict()
        for m in matches:
            gene2card[m['q_g']] = m['card']
        return gene2card


    def get_mini_contig(gene_name, window=3):

        # finds the contig containing the gene, and retreives the whole contig array
        match = col_neighs.find_one(
            {"genes.g": gene_name},
            {"c":1, "genes":1})

        if match:

            # Fix unordered contig problem > not needed anymore
            sorted_genes = sorted(match['genes'], key=lambda x: x['s'])
            for pos, g in enumerate(sorted_genes):
                g['p'] = pos

            # extract region from the whole contig
            anchor = next(pos for pos, g in enumerate(sorted_genes) if g['g'] == gene_name)
            start = max(0, anchor-window)
            end = anchor+window
            return sorted_genes[start:end+1]
        else:
            return []


    def get_cards(names):
        matches = col_cards.find({'q_g': {'$in': names} })
        gene2card = defaultdict()
        for m in matches:
            gene2card[m['q_g']] = m['card']
        return gene2card


    def add_to_cons_score(rel_pos,db_name,annotation,n_strand,anchor_strand,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict,i,distance):
        cons_per_position[rel_pos][db_name][annotation] += 1
        if n_strand != anchor_strand:
            strand_flag_per_pos[rel_pos][annotation] += 1
            neigh_number_per_pos_neg[rel_pos][annotation].append(str(i))
        neigh_number_per_pos[rel_pos][annotation].append(str(i))
        distances_dict[rel_pos][annotation].append(str(distance))
        return cons_per_position,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict


    fam,n, members = map(str.strip, line.split('\t'))
    nseqs = len(list(members.split(',')))

    # process each member of the family
    number_annots = defaultdict(lambda:defaultdict(lambda: Counter()))
    cog_descs = {}
    cons_per_position = defaultdict(lambda:defaultdict(lambda:Counter()))
    strand_flag_per_pos = defaultdict(lambda:Counter())
    neigh_number_per_pos = defaultdict(lambda:defaultdict(lambda:[]))
    neigh_number_per_pos_neg = defaultdict(lambda:defaultdict(lambda:[]))
    number_genes_per_pos = Counter()
    distances_dict = defaultdict(lambda:defaultdict(lambda:[]))

    for i,gene_entry in enumerate(members.split(',')):

        src,genome,gene,t = gene_entry.split('@')

        # First, give me neighbours and their positions/strands. The result includes the anchor
        mini_contig = get_mini_contig(gene, window=3)
        if mini_contig == []:
            continue

        # extract gene names from the mini contig
        mini_contig_genes = list([n['g'] for n in mini_contig])

        # query their annotations
        gene2annot = get_emapper_annotations(mini_contig_genes)
        #gene2card = get_cards(mini_contig_genes)


        # get strand of query gene and turn around array if -
        anchor_strand = [n['o'] for n in mini_contig if n['g'] == gene]
        if anchor_strand[0] == '-':
             mini_contig_genes = mini_contig_genes[::-1]

        # get position of query genes to have as reference
        pos_anchor = mini_contig_genes.index(gene)

        # get distances between genes
        distances = get_distances (mini_contig,pos_anchor,mini_contig_genes,anchor_strand)

        #print (fam,mini_contig_genes)
        for pos,neighbour_gene in enumerate(mini_contig_genes):
            rel_pos = pos - pos_anchor
            number_genes_per_pos[rel_pos] += 1
            n_strand = [n['o'] for n in mini_contig if n['g'] == neighbour_gene]

            distance = 0
            if rel_pos < 0:
                distance = distances[rel_pos][rel_pos+1]
            elif rel_pos >0:
                distance = distances[rel_pos-1][rel_pos]

            if neighbour_gene in gene2card:
                annotation = gene2card[neighbour_gene]
                db_name = "CARD"
                cons_per_position,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict = add_to_cons_score(rel_pos,db_name,annotation,n_strand,anchor_strand,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict,i,distance)

            if neighbour_gene in gene2annot:
                for db_name in gene2annot[neighbour_gene]:
                    if db_name in ['ogs','kpath','pfam','kos','cazy','kmods']:
                        for annotation in gene2annot[neighbour_gene][db_name]:
                            cons_per_position,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict = add_to_cons_score(rel_pos,db_name,annotation,n_strand,anchor_strand,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict,i,distance)

                    elif db_name == 'pname':
                        annotation = gene2annot[neighbour_gene][db_name]
                        cons_per_position,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict = add_to_cons_score(rel_pos,db_name,annotation,n_strand,anchor_strand,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict,i,distance)

            else: # also include genes wo annotation
                db_name = 'unknown'
                annotation = 'unknown'
                cons_per_position,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict = add_to_cons_score(rel_pos,db_name,annotation,n_strand,anchor_strand,strand_flag_per_pos,neigh_number_per_pos_neg,neigh_number_per_pos,distances_dict,i,distance)

    array_to_report = []
    for pos in cons_per_position:
        for db in cons_per_position[pos]:
            for annot,number in cons_per_position[pos][db].items():
                array_to_report.append ('\t'.join([fam,str(pos),db,annot,str(number),str(number_genes_per_pos[pos]),str(nseqs),','.join(neigh_number_per_pos[pos][annot]),','.join(neigh_number_per_pos_neg[pos][annot]),','.join(distances_dict[pos][annot]),str(strand_flag_per_pos[pos][annot])]))

    return array_to_report



for ln, line in enumerate(open(sys.argv[1])):
    print('\n'.join(process_fam(line)))
