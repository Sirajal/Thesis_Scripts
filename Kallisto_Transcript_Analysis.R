# source("http://bioconductor.org/biocLite.R")

#For filtering
library(tidyverse)
#differential transcript
library(sleuth)
#annotation
library(biomaRt)

#Sleuth does not do contrast like DESeq2 so need to only read heart or embryo
#samples are read in. That means the first bit is what tissue to run
# tissue_in = "Heart"
tissue_in = "Embryo"

#This is the same as the gene level for making sample table
#extra part for filtering tissue
pheno_in = list.files("../Data/Kallisto/", full.names = T) %>% # List the salmon output folders
  data.frame(Dir = .) %>% # Make into a data frame
  mutate(Tissue = ifelse(grepl("heart", Dir), "Heart", "Embryo"), # Make a column that gets the type of tissue
  Type = ifelse(grepl("WT", Dir), "WT", "Mutant"), # Make a column that lists the WT or Mutant samples
  Rep = rep(1:3, 4), # Make a replicate number
  Var = paste0(Tissue, "_", Type), #
  sample = paste0(Var, Rep)) %>%
  dplyr::select(sample, Var, Tissue, Type, path = Dir) %>%
  mutate_each(funs(as.vector), path) %>%
  filter(Tissue == tissue_in)


#make sleuth object as manual shows
so = sleuth_prep(pheno_in)
#make full model as manual shows
so = sleuth_fit(so, ~Var, 'full')
#make reduced model as manual shows
so = sleuth_fit(so, ~1, 'reduced')
#run test as manual shows
so = sleuth_lrt(so, 'reduced', 'full')
#get results as manual shows
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = F)
#filter results for adjusted p val less than 0.05
sleuth_significant = filter(sleuth_table, qval < 0.05)


#This works in the same way as the gene level
#difference that the filter is by transcript ID instead of gene id
ensembl = useMart("ensembl", dataset = "drerio_gene_ensembl")
anno = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","zfin_id_symbol","hgnc_symbol","description","entrezgene"),
        filters = "ensembl_transcript_id", values = sleuth_significant$target_id, mart = ensembl)


#add annotation to result
sleuth_significant <- sleuth_significant %>% left_join(anno, by = c("target_id"="ensembl_transcript_id"))
# Write result to csv file
write_csv(sleuth_significant, path = paste0("../Results/Kallisto_Transcript/TranscriptDE_", tissue_in, "_SignificantHits.csv"))
# graph jupa
plot_bootstrap(so, "ENSDART00000146512", units = "est_counts", color_by = "Var") + theme_bw()
