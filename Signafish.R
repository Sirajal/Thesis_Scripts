source("http://bioconductor.org/biocLite.R")
library(tidyverse)
library(biomaRt)

#signafish database from website http://signafish.org/download
signafish_db = as.data.frame(read_delim("02-02-2016-signafish-1.0.0-ZKqlql.csv", delim = ";"))
#get pathways in database
pathways = unique(na.omit(signafish_db$target_pathways))
pathways = strsplit(as.vector(pathways),',')
pathways = unique(unlist(pathways))

#get result
#Genes_In = read_csv("../Results/Kallisto_Gene/DE_GeneLevel_Kallisto_Results_Embryo_NoFC.csv")
#Genes_In = read_csv("../Results/Kallisto_Gene/DE_GeneLevel_Kallisto_Results_Heart_NoFC.csv")
#Genes_In = read_csv("../Results/Salmon_Gene/DE_GeneLevel_Salmon_Results_Embryo_NoFC.csv")
#Genes_In = read_csv("../Results/Salmon_Gene/DE_GeneLevel_Salmon_Results_Heart_NoFC.csv")
#Genes_In = read_csv("../Results/Stringtie_Gene/DE_GeneLevel_Stringtie_Results_Embryo_NoFC.csv")
Genes_In = read_csv("../Results/Stringtie_Gene/DE_GeneLevel_Stringtie_Results_Heart_NoFC.csv")


results = c()
#loop pathways and run hypergeometric test
for(i in pathways) {
  #database for pathway
  pathwaygenes = signafish_db[signafish_db$target_pathways == i,]
  #genes in signafish pathway
  wballs = unique(na.omit(pathwaygenes$target_speciesID))
  #genes in signafish not in pathway
  bballs = setdiff(unique(na.omit(pathwaygenes$target_speciesID)),wballs)
  #number of gene in pathway and differential expression
  qin = length(intersect(Genes_In$ensembl_gene_id, wballs))
  #run test
  hypergeom = phyper(q=qin,m=length(wballs),n=length(bballs),k=length(Genes_In$ensembl_gene_id),lower.tail = TRUE,log.p = FALSE)
  #add to results
  results = rbind(results, data.frame(pathway = i, pval = (1-hypergeom)))
}

#adjust p val
results$PAdj = p.adjust(results$pval, method = "BH")
#write results to excel file
write_csv(results, path = "../Results/Signafish/Stringtie/Signafish_Heart_Stringtie.csv")
