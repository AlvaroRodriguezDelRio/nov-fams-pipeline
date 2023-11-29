library(data.table)
library(caret)
library(dplyr)
library(coin)
library(gtools)
library(tidyr)


# load abundance data per novel family
data_n = read.csv("Abs_per_family.csv", header = T, stringsAsFactors=F)

# load metadata
samples = read.table("metadata.tsv",header = T)

# include colonoscopy into study for correcting in pvalue function for calculating p-values
samples <- samples %>%
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(Study!='CN.CRC', Study,
                      paste0(Study, '_', Sampling_rel_to_colonoscopy)))

samples$Study = samples$block

# merge abundance data and metadata
data_m = merge(data, samples, by = "sample")

#  Wilcoxon test blocked by population and colonoscopy for each novel family
mean_comparison = function(x) p.adjust(pvalue(wilcox_test(x ~ Group | Study, data=data_m)),method = 'fdr',n = nfams)
p_values = as.matrix(lapply(data_m[2:(ncol(data_m)-2)],mean_comparison))
