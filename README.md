# nov-fams-pipeline

Here we present the pipeline for computing novel gene families from the proteomes of a collection of genomes, and how to calculate their genomic context conservation. 

The scripts presented here assume that the gene names are formatted as >genome_source_of_isolation@genome_name@gene_name@domain|phylum 

## Deep homology-based protein clustering

Run mmseqs for calculating the gene families on the proteomes of interest (we used the  ```--min-seq-id 0.3 -c 0.5 --cov-mode 1 --cluster-mode 2 -e 0.001``` parameter combination).

## Detection of protein clusters specific from uncultivated taxa

Map gene families against reference databases for isolating those exclusive on uncultivated taxa.

- EggNOG with eggnog-mapper. We used the following parameter combination: ```emapper.py -m diamond --itype proteins --no_file_comments --cpu 5 -i multifasta.faa -o mappings.tab```). Eggnog annotations will also be useful for calculating the genomic context conservation of the gene families. We considered as significant any hit with E-value < 1e-3.

- PfamA with the command ```hmmsearch --cpu 10 --tblout mappings.tab /data/Pfam/Pfam-32-A/Pfam-A.hmm multifasta.faa```. We considered as significant any hit with E-value < 1e-5

- PfamB: with the command ```/scratch/alvaro/DEEM/analysis/mappings_db/pfamB/scripts/hmmsearch.sh```. We considered as significant any hit with E-value < 1e-5.

- RefSeq with the command ```diamond blastx -d /data/RefSeq/refseq.dmnd -q multifasta.cds -o results.tab --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp qlen slen --sensitive```. Hits with an E-value < 1e-3 and query coverage > 50% were considered significant

We considered as novel gene families those with no significant hits to any of these databases.

## Reconstructing the genomic context of novel families 

For reconstructing genomic contexts: 

- 

get a list with all the genes within each genome and contig (python /scratch/alvaro/DEEM/analysis/neighs/build_neighbours_per_contig/scripts/neighs_per_contig.py paths_gffs.txt > neighs_per_contig.tab  (input is, for instance,  /scratch/alvaro/DEEM/analysis/data/paths_gff.txt)).
update them to a mongo collection: python /scratch/alvaro/DEEM/analysis/build_db/scripts/neigh2json.py neighs_per_contig.tab [...] | mongoimport --host fat01 -d XXX -c neighs --drop. If you want to reconstruct the genomic context of the extended families, you also need to input the neigh_per_contig.tab file for the 169k genomes (/data/jhc/cold/MAGs/GEM-fixed_genomes/neighs_per_contig.tab /data/jhc/cold/MAGs/GMGC-v1/high-q/neighs_per_contig.r.tab /data/jhc/freezer/public/GTDB/GTDB_rev95/neighs_per_contig.no_dups.tab /data/jhc/cold/MAGs/Ocean-v1/neighs_per_contig.tab /data/jhc/cold/MAGs/UHGG-v1/neighs_per_contig.tab). 
Get name of all neighbors of the genes in the novel families (only needed if reconstructing the genomic context after expanding, for getting the emapper annotations of the genes from the 169k genomes) with python  /scratch/alvaro/DEEM/analysis/neighs/scripts/get_names_neigh_genes.py > neigh_genes.txt (gets as input gene family composition after extending). 
Get mongo collection with the functional annotation of all neighbors (including those in the 169k genomes, this script looks for the annotations of the genes in the neigh_genes.txt list in a database containing the annotations of those genomes. if given an empty list, it will only get the funct annotation from the eggnogmapper_out.tab file): python /scratch/alvaro/DEEM/analysis/build_db/scripts/emapper.py eggnogmapper_out.tab | mongoimport --host fat01 -d DATABASE_NAME -c emapper2 --drop (emapper-2.1.5 output). 
get genomic context conservation:
precompute: python  /scratch/alvaro/DEEM/analysis/neighs/scripts/score_per_pos.strand.unknown_genes.space.py extended_gene_family_composition.tab > scores.tab; 
json: python /scratch/alvaro/DEEM/analysis/neighs/scripts/score_per_pos.strand2json.espace.py scores.tab > scores.json;
tabular format: python/scratch/alvaro/DEEM/analysis/neighs/scripts/score_per_pos.strand2table.py (fields: family name, db, functional_term, position, score, % cont strand, % contrary strand in between the novel genes and the genes with the functional term, % genes separated more than 100nts in between the novel gene and the neighbors, description).
 The json needs to be uploaded to the mongo database for visualizing the genomic context: cat conserved_context.json | mongoimport --host fat01 -d XXXXXX -c og_neigh_scores --drop

Taxonomic coverage and specificity: python /scratch/alvaro/DEEM/analysis/taxonomic_analyses/scripts/get_coverage_sp.raw.py family_composition.tab > sp_cov_per_fam_per_lin.tab. The script also reads from tax annotation per genome files (needs to be changed in the script for your custom genomes, the tax annota for the 169k genomes file is already read by the script. 
Get individual fasta per family, alignments and trees
file per fam (subdirs need to be created beforehand): python  /scratch/alvaro/DEEM/analysis/file_per_extended_family/scripts/get_individual_fas_per_fam.py multifasta.faa [/data/jhc/cold/MAGs/combined/microbial_genomes-v1.proteins.dups_reformatted.faa if want to incorporate genes from extended families]
alignments and trees: take file of fasta files: sbatch /scratch/alvaro/DEEM/analysis/file_per_extended_family/scripts/alignments.sh
get tree jsons: python /scratch/alvaro/DEEM/analysis/build_db/scripts/get_tree_json.py paths_trees.txt | mongoimport --host fat01 -d DATABASE_name -c trees --drop
Signal peptides / TM domains: 
calculate (all takes as input paths with individual fasta per fam): signalP: 
/scratch/alvaro/DEEM/analysis/sequence_signatures/signalP/scripts/signalP.gram_neg.sh and /scratch/alvaro/DEEM/analysis/sequence_signatures/signalP/scripts/signalP.gram_pos.sh
tmhmm:   /scratch/alvaro/DEEM/analysis/sequence_signatures/tmhmm/scripts/tmhmm.sh (denbi)
process:  /scratch/alvaro/DEEM/analysis/sequence_signatures/tmhmm/scripts/tmhmm2json.py (denbi) and /scratch/alvaro/DEEM/analysis/sequence_signatures/signalP/scripts/signalP2json.py (cbgp). This generates jsons which can be uploaded to the mongo database (cat signalp.json | mongoimport --host fat01 -d DATABASSE_NAME -c signalp --drop, cat tmhmm.json| mongoimport --host fat01 -d DATABASE_NAME -c tm --drop)
For visualizing Families in the online resource by Jorge, the database also need to have the following collections: 
fam2info, with python/scratch/alvaro/DEEM/analysis/build_db/scripts/fam2info.noemapper.py |   mongoimport --host fat01 -d DATABASE_NAME -c faminfo - -drop
fam2members: python /scratch/alvaro/DEEM/analysis/build_db/scripts/fam2members2json.py |  mongoimport --host fat01 -d DATABASE_NAME -c fam2members --drop
genome 2 taxonomy: python /scratch/alvaro/DEEM/analysis/build_db/scripts/genometax2json.py |  mongoimport --host fat01 -d DATABASE_NAME -c genome_taxonomy --drop (change genome2taxonomy.names.tab to /data/jhc/cold/MAGs/combined/genome2taxonomy.r202.tab). 
