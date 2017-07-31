
source("http://bioconductor.org/biocLite.R")
library(tidyverse)
library(biomaRt)
#to plot https://github.com/hms-dbmi/UpSetR
library(UpSetR)


#get result of differential expression of genes
salmon.deseq2.results.embryo = read_csv("../Results/Salmon_Gene/DE_GeneLevel_Salmon_Results_Embryo_FC1_5.csv")
kallisto.deseq2.results.embryo = read_csv("../Results/Kallisto_Gene/DE_GeneLevel_Kallisto_Results_Embryo_FC1_5.csv")
stringtie.deseq2.results.embryo = read_csv("../Results/Stringtie_Gene/DE_GeneLevel_Stringtie_Results_Embryo_FC1_5.csv")
salmon.deseq2.results.heart = read_csv("../Results/Salmon_Gene/DE_GeneLevel_Salmon_Results_Heart_FC1_5.csv")
kallisto.deseq2.results.heart = read_csv("../Results/Kallisto_Gene/DE_GeneLevel_Kallisto_Results_Heart_FC1_5.csv")
stringtie.deseq2.results.heart = read_csv("../Results/Stringtie_Gene/DE_GeneLevel_Stringtie_Results_Heart_FC1_5.csv")
results.list = list()
results.list[["salmon.deseq2.results.embryo"]] = salmon.deseq2.results.embryo
results.list[["kallisto.deseq2.results.embryo"]] = kallisto.deseq2.results.embryo
results.list[["stringtie.deseq2.results.embryo"]] = stringtie.deseq2.results.embryo
results.list[["salmon.deseq2.results.heart"]] = salmon.deseq2.results.heart
results.list[["kallisto.deseq2.results.heart"]] = kallisto.deseq2.results.heart
results.list[["stringtie.deseq2.results.heart"]] = stringtie.deseq2.results.heart

#make dataframe for plot. column is differential expression file
#row is gene. 1 for differential express 0 for not
upsetdf = data.frame(ensembl_gene_id = "")
#loop differential expression results
#make dataframe of 1 and ensembl_gene_id
#full join keeps missing gene not in another list
for(i in 1:length(results.list)) {
  tmp = data.frame(ensembl_gene_id = results.list[[i]]$ensembl_gene_id, DE = 1)
  colnames(tmp) = c("ensembl_gene_id", names(results.list)[i])
  upsetdf = full_join(upsetdf, tmp)
}
#delete row 1 all 0
upsetdf = upsetdf[-1,]
#missing gene is na, make 0
upsetdf[is.na(upsetdf)] = 0
#make upset plot
upset(upsetdf, sets = colnames(upsetdf)[-1])



#get result of differential expression of transcript
salmon.sleuth.results.embryo = read_csv("../Results/Salmon_Transcript/TranscriptDE_Embryo_SignificantHits.csv")
kallisto.sleuth.results.embryo = read_csv("../Results/Kallisto_Transcript/TranscriptDE_Embryo_SignificantHits.csv")
salmon.sleuth.results.heart = read_csv("../Results/Salmon_Transcript/TranscriptDE_Heart_SignificantHits.csv")
kallisto.sleuth.results.heart = read_csv("../Results/Kallisto_Transcript/TranscriptDE_Heart_SignificantHits.csv")
results.list = list()
results.list[["salmon.sleuth.results.embryo"]] = salmon.sleuth.results.embryo
results.list[["kallisto.sleuth.results.embryo"]] = kallisto.sleuth.results.embryo
results.list[["salmon.sleuth.results.heart"]] = salmon.sleuth.results.heart
results.list[["kallisto.sleuth.results.heart"]] = kallisto.sleuth.results.heart

#make dataframe for plot. column is differential expression file
#row is gene. 1 for differential express 0 for not
upsetdf = data.frame(ensembl_gene_id = "")
#loop differential expression results
#make dataframe of 1 and ensembl_gene_id
#full join keeps missing gene not in another list
for(i in 1:length(results.list)) {
  tmp = data.frame(ensembl_gene_id = results.list[[i]]$ensembl_gene_id, DE = 1)
  colnames(tmp) = c("ensembl_gene_id", names(results.list)[i])
  upsetdf = full_join(upsetdf, tmp)
}
#delete row 1 all 0
upsetdf = upsetdf[-1,]
#missing gene is na, make 0
upsetdf[is.na(upsetdf)] = 0
#make upset plot
upset(upsetdf, sets = colnames(upsetdf)[-1])



#get result of differential expression of exons
#make unique id for plot exon number and ensembl_gene_id
dexseq.embryo = read_csv("../Results/DEXSeq_Embryo/dexseq_results.csv") %>% mutate(uniqueid = paste0(groupID, featureID))
dexseq.heart = read_csv("../Results/DEXSeq_Heart/dexseq_results.csv") %>% mutate(uniqueid = paste0(groupID, featureID))
results.list = list()
results.list[["dexseq.heart"]] = dexseq.heart
results.list[["dexseq.embryo"]] = dexseq.embryo

#make dataframe for plot. column is differential expression file
#row is gene. 1 for differential express 0 for not
upsetdf = data.frame(uniqueid = "")
#loop differential expression results
#make dataframe of 1 and ensembl_gene_id
#full join keeps missing gene not in another list
for(i in 1:length(uniqueid)) {
  tmp = data.frame(ensembl_gene_id = results.list[[i]]$uniqueid, DE = 1)
  colnames(tmp) = c("uniqueid", names(results.list)[i])
  upsetdf = full_join(upsetdf, tmp)
}
#delete row 1 all 0
upsetdf = upsetdf[-1,]
#missing gene is na, make 0
upsetdf[is.na(upsetdf)] = 0
#make upset plot
upset(upsetdf, sets = colnames(upsetdf)[-1])
