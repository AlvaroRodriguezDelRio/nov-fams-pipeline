# Pipeline for novel family computation and analysis

Here we present the pipeline for computing novel gene families from the proteomes of a collection of genomes, and how to calculate their genomic context conservation.

![Pipeline for dealineating novel gene families exclusive of uncultivated taxa](Capture.PNG)

The scripts presented here assume that the gene names are formatted as ```>genome_source_of_isolation@genome_name@gene_name@domain|phylum```.

## Deep homology-based protein clustering

Run mmseqs for calculating the gene families on the concatenated proteomes of the genomes of interest (we used the  ```--min-seq-id 0.3 -c 0.5 --cov-mode 1 --cluster-mode 2 -e 0.001``` parameter combination).

## Detection of protein clusters specific from uncultivated taxa

Map gene families against reference databases for isolating those exclusive on uncultivated taxa.

- EggNOG with eggnog-mapper. We used the following parameter combination: ```emapper.py -m diamond --itype proteins --no_file_comments -i multifasta.faa -o mappings.tab```. We considered as significant any hit with E-value < 1e-3. Eggnog annotations will also be useful for calculating the genomic context conservation of the gene families.

- PfamA with the command ```hmmsearch --cpu 10 --tblout mappings.tab /data/Pfam/Pfam-32-A/Pfam-A.hmm multifasta.faa```. We considered as significant any hit with E-value < 1e-5. 

- PfamB: with the command ```/scratch/alvaro/DEEM/analysis/mappings_db/pfamB/scripts/hmmsearch.sh```. We considered as significant any hit with E-value < 1e-5. !!!!!!!!!!!!!!!!!!!!

- RefSeq with the command ```diamond blastx -d /data/RefSeq/refseq.dmnd -q multifasta.cds -o results.tab --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp qlen slen --sensitive```. Hits with an E-value < 1e-3 and query coverage > 50% were considered significant

We considered as novel gene families those with no significant hits to any of these databases.

## Filtering novel gene families 

For dealineating the novel gene family predictions of higher quality, we conducted a series of analysis, for which we needed to create individual fasta files with the i) CDS sequences and ii) Protein sequences of each gene family. 

- For calculating conserved domains, multiple sequence alignments need to be computed on each family. For such purpose, we used clustalO: 

``` clustalo command```

After having collected all the protein alignments, domain conservation can be calculated by running:

 ```python XXXXXX ```
 
- For discarding viral sequences from the novel gene familes, we mapped the protein sequences against the PVOGs database (https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/home.html). 

``` COMMAND ```

We considered hits with E-value < 1e-5 and minimum coverage of 50% as significant, and discarded families with significant hits.

- For discarding sporious sequences, the protein sequences can be mapped against the Antifam database (https://ftp.sanger.ac.uk/pub/databases/Pfam/AntiFam/). 

``` COMMAND```

We discarded families with any hit with E-value < 1e-5. 

- For calculating the dN/dS and coding probability of each gene family, we built protein alignments, back translated them to nucleotides, and calculated gene family trees with the ETE toolkit (http://etetoolkit.org/): 

```ETE command```

Later, we ran hyphy BUSTED (http://vision.hyphy.org/) for calculating the dN/dS of each gene family:

```BUSTED command```

And RNAcode (https://github.com/ViennaRNA/RNAcode) for calculating their coding probability:

```RNAcode command```

## Reconstructing the genomic context of novel gene families 

For reconstructing genomic contexts: 

- Get an ordered list of the genes within each contig with ```python neighs_per_contig.py paths_gffs.txt > neighs_per_contig.tab```

-----(input is, for instance,  /scratch/alvaro/DEEM/analysis/data/paths_gff.txt)).
------ /scratch/alvaro/DEEM/analysis/neighs/build_neighbours_per_contig/scripts/neighs_per_contig.py

- Update this data to a Mongo collection with  ```python /scratch/alvaro/DEEM/analysis/build_db/scripts/neigh2json.py neighs_per_contig.tab | mongoimport --host fat01 -d XXX -c neighs --drop```. CHANGE FOR ALSO LOADING DATA INTO MEMORY????

- Create a Mongo collection with the functional annotation of the neighbors of the members of each gene family: ```python /scratch/alvaro/DEEM/analysis/build_db/scripts/emapper.py eggnogmapper_out.tab | mongoimport --host fat01 -d DATABASE_NAME -c emapper2 --drop (emapper-2.1.5 output).```

- Calculate genomic context conservation:

1) Precompute:
```python  /scratch/alvaro/DEEM/analysis/neighs/scripts/score_per_pos.strand.unknown_genes.space.py extended_gene_family_composition.tab > scores.tab``` 

2) Genomic context conservation in json format.
```python /scratch/alvaro/DEEM/analysis/neighs/scripts/score_per_pos.strand2json.espace.py scores.tab > final_scores.json```

3) Genomic context conservation in tabular format
```python/scratch/alvaro/DEEM/analysis/neighs/scripts/score_per_pos.strand2table.py > final_scores.tab```. The fields in the table are: family name, db, functional_term, position, score, % cont strand, % contrary strand in between the novel genes and the genes with the functional term, % genes separated more than 100nts in between the novel gene and the neighbors, description.

## Gene family taxonomic coverage and specificity

For calculating the taxonomic coverage and specificity for each gene family on each taxonomic group, use:

```python /scratch/alvaro/DEEM/analysis/taxonomic_analyses/scripts/get_coverage_sp.raw.py family_composition.tab > sp_cov_per_fam_per_lin.tab.```

The script also reads from taxonomic annotation per genome files (needs to be changed in the script for your custom genomes, the tax annota for the 169k genomes file is already read by the script. #CHANGE IN SCRIPT FOR READING STDIN#

