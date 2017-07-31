#Load Libraries
library(tidyverse)
library(DEXSeq)

#Set working directory
setwd("~/Documents/RNA_Seq_Project/Scripts/")
#sample type
sample.in = "Embryo"
#Get file names and make data frame of samples
pheno.in = list.files("../Data/DEXSeq/", full.names = T) %>% # List the DEXSeq output files
  data.frame(File = .) %>% # Make into a data frame
  mutate(Tissue = ifelse(grepl("heart", File), "Heart", "Embryo"), # Make a column that gets the type of tissue
         Type = ifelse(grepl("WT", File), "WT", "Mutant"), # Make a column that lists the WT or Mutant samples
         Var = paste0(Tissue, "_", Type)) %>%
  arrange(Var) %>% # Arrange by Var so that they're in the correct order for replicate numbers
  mutate(Rep = rep(1:3, 4), # Make a replicate number)
         sample = paste0(Var, Rep)) %>%
  dplyr::select(sample, Var, Tissue, Type, path = File) %>%
  mutate_each(funs(as.vector), path) %>%
  filter(Tissue == sample.in) #Only analyse the Embryo samples

# Make a DEXseq dataset as the manual shows
dxd = DEXSeqDataSetFromHTSeq(countfiles = pheno_in$path,
                             sampleData = pheno_in,
                             design = ~ sample + exon + Type:exon,
                             flattenedfile = "../Reference_Files/Danio_rerio.GRCz10.88_Chr_DEXSeq.gff" )
# get size factors as manual shows
dxd = estimateSizeFactors(dxd)
# model fit as manual shows
dxd = estimateDispersions(dxd)
# run differential expression test as manual shows
dxd = testForDEU(dxd)
# calculate fold changes
dxd = estimateExonFoldChanges( dxd, fitExpToVar="Type") #Run this next
# create results
dxr = DEXSeqResults(dxd)
# plotDEXSeq(dxr, "", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
# save output as RData file because long time to calculate
save(dxd, dxr, file = paste0("../Preprocessed_Data/DEXSeq_", sample.in, ".RData"), compress = T)
# make html output plots as manual shows
# DEXSeqHTML(dxr, FDR=0.05, color=c("#FF000080", "#0000FF80"), file = paste0("../Results/DEXSeq_", sample.in, ".html"))



# load(paste0("../Preprocessed_Data/DEXSeq_", sample.in, ".RData"))
# make html output plots as manual shows
DEXSeqHTML(dxr, FDR=0.05, color=c("#FF000080", "#0000FF80"),
           fitExpToVar="Type", path = "./",
           file = paste0("./DEXSeq_", sample.in, ".html"))

# get significant exons adjusted p value less than 0.05, and greater than 1.5 fold change
dexseq_significant = dxr %>% as.data.frame %>%
  filter(padj < 0.05, abs(log2fold_WT_Mutant) > log2(1.5)) %>%
  arrange(desc(abs(log2fold_WT_Mutant))) %>%
  mutate(ensembl_gene_id = substr(groupID, 1, 18))

# get ensembl annotation using biomart, shown in manual
ensembl = useEnsembl("ensembl", dataset = "drerio_gene_ensembl")
anno = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "zfin_id_symbol"),
                 filters = "ensembl_gene_id", values = dexseq_significant$ensembl_gene_id,
                 mart = ensembl)

# get significant exons and annotate
dexseq_significant = dexseq_significant %>% left_join(anno) %>%
  dplyr::select(groupID,featureID,log2fold_WT_Mutant,padj,ensembl_gene_id,
                zfin_id_symbol, hgnc_symbol, description)
# write results to file
write_csv(dexseq_significant, path = "../Results/DEXSeq_Embryo/dexseq_results.csv")
