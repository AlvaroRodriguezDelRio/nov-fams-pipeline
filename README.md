# Pipeline for functionally and evolutionarily significant novel (FESNov) gene families computation and characterization

Here we present the pipeline for computing novel gene families from the proteomes of a collection of genomes, and how to calculate their genomic context conservation.

![Pipeline for dealineating novel gene families exclusive of uncultivated taxa](Pipeline.png)

The starting point of this pipeline consist a the concatenated fasta file with the gene predictions for the genomes of interest. The scripts presented here assume that the gene names are formatted as ```>genome_source_of_isolation@genome_name@gene_name@domain|phylum```.

## Requirements

### Software

- MMseqs2: https://github.com/soedinglab/MMseqs2
- diamond: https://github.com/bbuchfink/diamond
- HMMER: http://hmmer.org/
- clustalO: http://www.clustal.org/omega/
- FastTree: http://www.microbesonline.org/fasttree/
- ETE3: http://etetoolkit.org/download/
- HyPhy: https://www.hyphy.org/
- RNAcode: https://github.com/ViennaRNA/RNAcode
- mongoDB: https://www.mongodb.com/
- GTDB-Tk: https://github.com/Ecogenomics/GTDBTk

### Databases

- RefSeq: https://www.ncbi.nlm.nih.gov/refseq/
- EggNOG: http://eggnog5.embl.de/#/app/home
- Pfam: http://pfam.xfam.org/
- PVOGs: https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/home.html

## Deep homology-based protein clustering

Run MMseqs2 for calculating the gene families on the concatenated proteomes of the genomes of interest (we used the  ```--min-seq-id 0.3 -c 0.5 --cov-mode 1 --cluster-mode 2 -e 0.001``` parameter combination).

## Detection of protein clusters specific from uncultivated taxa

We then mapped the gene families against the following reference databases for isolating those exclusive on uncultivated taxa.

- EggNOG with eggnog-mapper. We used the following parameter combination: ```emapper.py -m diamond --itype proteins --no_file_comments -i multifasta.faa -o mappings.tab```. We considered as significant any hit with E-value < 1e-3. 

- PfamA with the command ```hmmsearch --tblout mappings.tab /data/Pfam/Pfam-32-A/Pfam-A.hmm multifasta.faa```. We considered as significant any hit with E-value < 1e-5. 

- All the hmms within the PfamB collection with the command ```hmmsearch --tblout mappings.tab PfamB.hmm multifasta.faa```. We considered as significant any hit with E-value < 1e-5. 

- RefSeq with the command ```diamond blastx -d /data/RefSeq/refseq.dmnd -q multifasta.cds -o results.tab --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp qlen slen --sensitive```. Hits with an E-value < 1e-3 and query coverage > 50% were considered significant

We considered as gene families exclusive from uncultivated taxa those with no member with significant hits to any of these databases.

## Unknown gene familiy from uncultivated taxa filtering 

For dealineating the unknown gene family predictions of higher quality, we conducted a series of filtering steps, for which we needed to create individual fasta files with the i) CDS sequences and ii) Protein sequences of each gene family. 

- For calculating conserved domains, multiple sequence alignments need to be computed on each family. For such purpose, we used clustal omega (http://www.clustal.org/omega/): 

```clustalo -i gene_family_protein_fasta.faa -o gene_family_protein_fasta.alg.faa```

After having collected all the protein alignments, domain conservation can be calculated by running:

 ```python calculation_conserved_domain.py paths_algs.txt```

 We discarded gene families with a conserved domain shorter than 20 residues.
 
- For discarding viral sequences from the novel gene familes, we mapped the protein sequences against all the hmms within the PVOGs database (https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/home.html). 

```hmmsearch --tblout mappings.tab PVOG.hmm multifasta.faa```

We considered hits with E-value < 1e-5 and minimum coverage of 50% as significant, and discarded families with significant hits.

- For discarding sporious sequences, the protein sequences can be mapped against the Antifam database (https://ftp.sanger.ac.uk/pub/databases/Pfam/AntiFam/). 

```hmmsearch --cut_ga --tblout mappings_antifam.tab AntiFam.hmm  multifasta.faa```

We discarded families with any hit with E-value < 1e-5. 

- For calculating the dN/dS and coding probability of each gene family, we ran gene family trees on the protein alignments with FastTree (default options), and back translated the protein alignments to nucleotides with the ETE toolkit (http://etetoolkit.org/) 

```ete3 build -a gene_family_protein_fasta.faa -n gene_family_CDS_fasta.cds -o output_dir --nt-switch-threshold 0.0 --noimg -w clustalo_default-none-none-none```

Later, we ran 

-- hyphy BUSTED (http://vision.hyphy.org/) for calculating the dN/dS of each gene family:

```hyphy busted --alignment gene_family_CDS_fasta.alg.cds --tree gene_family.nw```

-- RNAcode (https://github.com/ViennaRNA/RNAcode) for calculating their coding probability. Before running RNAcode, alignments need to be changed to MAF format, which can be done in python by importing SeqIO from the Bio package (```SeqIO.parse(sys.argv[1], "fasta"); SeqIO.write(records, sys.argv[1]+'.maf', "maf")```)

```RNAcode gene_family_CDS_fasta.alg.cds.maf -t --stop-early -o RNAcode_out.tab```

We only considered in our analysis novel gene families with a conserved domain of at least 20 residues, with no significant homolgy in the pVOGs and Antifam databases, , with dN/dS < 0.5 and with RNAcode p-value < 0.05. 

## Reconstructing the genomic context of novel gene families 

For reconstructing the genomic context of the novel gene families, we followed the subsequent steps: 

- Get an ordered list of the genes within each contig with ```python genes_per_contig.py paths_gffs.txt > neighs_per_contig.tab```

The ```paths_gffs.txt``` file contains the gff files of all the genomes with the gene prediction information, for instance:

```
Chip-388_95C1R_METABAT_1       prokka  gene    85      783     .       +       .       ID=Chip-388_95C1R_METABAT_00001_gene;Name=tssA_1;gene=tssA_1;locus_tag=Chip-388_95C1R_METABAT_00001
Chip-388_95C1R_METABAT_1       prokka  mRNA    85      783     .       +       .       ID=Chip-388_95C1R_METABAT_00001_mRNA;Name=tssA_1;gene=tssA_1;locus_tag=Chip-388_95C1R_METABAT_00001
Chip-388_95C1R_METABAT_1       Prodigal:002006 CDS     85      783     .       +       0       ID=Chip-388_95C1R_METABAT_00001;Parent=Chip-388_95C1R_METABAT_00001_gene,Chip-388_95C1R_METABAT_00001_mRNA;eC_number=2.8.1.1;Name=tssA_1;db_xref=COG:COG2897;gene=tssA_1;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:D4GYM0;locus_tag=Chip-388_95C1R_METABAT_00001;product=Putative thiosulfate sulfurtransferase;protein_id=gnl|DEEM|Chip-388_95C1R_METABAT_00001
```

- Upload the ```neighs_per_contig.tab``` data to a Mongo (https://www.mongodb.com/) collection with  ```python neigh2json.py neighs_per_contig.tab | mongoimport --host HOST_NAME -d DATABASE_NAME -c neighs --drop```.  

- Create a Mongo collection with the functional annotation of the neighbors of the members of each gene family: ```python emapper2json.py eggnogmapper_out.tab | mongoimport --host HOST_NAME -d DATABASE_NAME -c emapper2 --drop (emapper-2.1.5 output)```.

  It is worth noting that creating Mongo databases may not be needed if working with a relatively low number of genomes. Then, the genomic context and gene functional information may be directly loaded into memory.  

- Calculate genomic context conservation:

1) Precompute:
```python neigh_cons_score.py gene_family_composition.tab > scores.tab```. The ```gene_family_composition.tab``` file is a tab delimeted file with 3 columns: the gene family name, number of members and coma-separated list of members. Then run ```python genomic_context_conservation_table.py scores.tab > final_scores.tab``` for generating a final report of genomic context conservation per gene family. The fields included in the final table are:

```
gebe_family_name
db (database from which the functional term comes from)
functional_term
relative_position (genomic position relative to the gene family members, e.g. -1 means gene immediate before the gene family members)
vertical score (% of neighboring genes in the position specified annnotated to the functional_term specified)
% of neighbors in the position & with the functional_term in contrary strand to the gene family members. For each position and functional term for which a score is provided, proportion of gene neighbors that are in contrary orientations than the gene family members. This informs of whether neighbors are in the same strand as the gene family members or not. 1 means that all neigbors in the position specified & showing the functional term are in contrary orientation. 
number of neighbors in contrary strand in between the novel genes and the genes in the position & with the functional_term. In case the position is not adjaccent to the gene family members (positions +-1), it may be the case that there are genes in between the gene family members and the genes in the position reported that are in contrary strand and break the operon structure. This field informs about this. 
% genes separated more than 100nts in between the novel gene family members and the genes in the position & with the functional_term, description of the functional term. This field represents the number of genes in between the target position and the gene family members which are more than 100nts away. Higher values indicate that separation between genes is generally higher than 100nts. 
```

This is how the table looks like:

```
GTDBiso@GB_GCA_003141455@PLMA01000055.1_15@d__Bacteria|p__Chloroflexota NOVPT00U        og      COG0627 -3      0.14285714285714285     0.0     0.0     0.0     Serine hydrolase involved in the detoxification of formaldehyde
GTDBiso@GB_GCA_003141455@PLMA01000055.1_15@d__Bacteria|p__Chloroflexota NOVPT00U        kpath   01130   2       0.7142857142857143      0.0     0.0     0.0
GTDBiso@GB_GCA_003141695@PLNE01000019.1_4@d__Bacteria|p__Dormibacterota NOV6H718        kpath   00020   -3      0.375   1.0     1.0      1.0     Citrate cycle (TCA cycle)
```

The ```genomic_context_confidence.py``` code was used to estimate the confidence of the KEGG pathway functional assignations based on genomic context, by measuring how different genomic architectures could correctly predict pathways on known function genes. The script reads from the ```final_scores.tab``` table, and also needs a table with confident KEGG pathway annotations per gene family (Fields: ```family_name,db from where the functional annotation come from,number of gene family members,funct_annot,n_genes_annot_as_funct_annotated_to_the_funct_annot_specicied,funct_annot_desc```). For each of these gene families, the script tests whether different genomic context conservation thresholds correcly predict the original KEGG pathway of the gene family.

## Gene family taxonomic coverage and specificity for synapomorphy discovery

For running this step, you need to compute the taxonomy of the genomes with GTDB-tk. For calculating the taxonomic coverage and specificity for each gene family on each taxonomic group, we used:

```python tax_cov_sp.py genome_tax_annotation.tab gene_family_composition.tab > sp_cov_per_fam_per_lin.tab.```. The ```genome_tax_annotation.tab``` file is a tab-delimeted file with 2 columns: genome name and their GTDB (https://gtdb.ecogenomic.org/) taxonomic annotations. ```sp_cov_per_fam_per_lin.tab``` is a tab-delimted file containing the following columns: gene family name, GTDB lineage, specificity, coverage, number of genomes within the GTDB lineage in the collection. 

## Sample discrimination and biomarker discovery 

For testing whether novel families could discriminate between conditions (i.e. control and colorectal cancer samples), we used the ```CRC prediction.R```, which reads an abundance table, with the abundance of each gene family in each sample, and calculates their prediction power using logistic models and machine learning algorithms.

For detecting novel gene families over or under represented in particular conditions, we ran the ```CRC_DA_analysis.r``` code, which reads a gene family abundance table per sample, and tests for differential abundance by calculating p-values by Wilcoxon test adjusted by the FDR method for multiple hypothesis testing. 
 
