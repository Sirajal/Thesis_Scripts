# source("https://bioconductor.org/biocLite.R") # Get bioconductor installer
# biocLite("tximport") # install a bioconductor

# read Kallisto
library(tximport)
# differential expression
library(DESeq2)
# read gtf (Gene to Transcript)
library(GenomicFeatures)
# For filtering
library(tidyverse)
library(reshape2)
# annotation
library(biomaRt)

#%>% pass the result
#paste0 function add word or variable together


#Reads in the kallisto file and make table to show samples
#first line list all output folder kallisto disrectory
#next line makes it in data frame
#next line starts the "mutate" part, which mean that we change column
#Tissue. Make col call Tissue, and if word "heart" is in file, call it "Heart", if it not, call "Embryo"
#Type. Same as tissue, but look for wt in file
#Rep. Replicate, number 1-3, four time. This is so individual filenames. deseq2 need individual
#Var. variable of interest in model. is Tissue, type column, with underline.
#sample. sample name, var and rep
#dplyr::select mean use select from tidyverse. order columns for deseq2
#mutate_each(funs(as.vector), path). path into vector
pheno_in = list.files("../Data/Kallisto/", full.names = T) %>% # List the kallisto output folders
data.frame(Dir = .) %>% # Make into a data frame
mutate(Tissue = ifelse(grepl("heart", Dir), "Heart", "Embryo"), # Make a column that gets the type of tissue
Type = ifelse(grepl("WT", Dir), "WT", "Mutant"), # Make a column that lists the WT or Mutant samples
Rep = rep(1:3, 4), # Make a replicate number
Var = paste0(Tissue, "_", Type), #
sample = paste0(Var, Rep)) %>%
dplyr::select(sample, Var, Tissue, Type, path = Dir) %>%
mutate_each(funs(as.vector), path)


#Take gtf file from ensembl and make it to a TxDb objects.
#This mean that we want the "GENEID" from the database
#columns in the right order
txdb = makeTxDbFromGFF("../Reference_Files/Danio_rerio.GRCz10.88.gtf", format = "gtf")
k = keys(txdb, keytype = "GENEID")
df = AnnotationDbi::select(txdb, keys = k , keytype = "GENEID" , columns = "TXNAME")
tx2gene = df[,c(2,1)]

#files_in make files for kallisto
#name files for kallisto sample name
#import kallisto files
files_in = paste0(pheno_in$path,"/abundance.h5")
names(files_in) = pheno_in$sample
txi = tximport(files_in, type = "kallisto" , tx2gene = tx2gene, varReduce = T)


#Make DESeq2Dataset for tximport
#Design is variable Var, as in manual
#Remove low counts, 60 reccommended in post. reduce noise
#Run DESeq model
#Normalised counts
dds = DESeq2::DESeqDataSetFromTximport(txi = txi, colData = pheno_in, design = ~ Var) # Make DESeq2 Dataset
dds_in = dds[rowSums(counts(dds)) > 60,]
dds_fit = DESeq(dds_in)
ncounts = counts(dds_fit, normalized = T)

#annotation biomart
#Get biomart code as in manual
#Get annotation attributes, filter on ensembl gene code
ensembl = useMart("ensembl", dataset = "drerio_gene_ensembl")
anno = getBM(attributes = c("ensembl_gene_id", "zfin_id_symbol","hgnc_symbol","description","entrezgene"),
        filters = "ensembl_gene_id", values = rownames(dds_fit), mart = ensembl)


#Principle Components Analysis (PCA)
#deseq2 method for making counts log
#deseq2 method for plotting PCA
rld = rlogTransformation(dds_fit)
plotPCA(rld, intgroup=c("Tissue", "Type"))


#Get results, like in manual
#get result for different FC - 0,1.5,2,3,4
#For loop
foldchange_parameters = c("NoFC" = 0, "FC1_5" = 1.5, "FC2" = 2, "FC3" = 3, "FC4" = 4)
for(i in 1:length(foldchange_parameters)) {
  #if FC is 0, make fc_in 0, as log2 0 not work
  if(foldchange_parameters[i] == 0) {
    fc_in = 0
  } else {
  #log2 FC
    fc_in = log2(foldchange_parameters[i])
  }
  #Name of FC
  fc_name = names(foldchange_parameters)[i]

  #get result of embryo mutant wt
  #add rowname for ensembl
  #filter FC and adjusted p val. abs mean no negative
  #add annotation
  res_first = as.data.frame(results(dds_fit, contrast = c("Var", "Embryo_Mutant", "Embryo_WT"))) %>%
      add_rownames("ensembl_gene_id") %>%
      filter(abs(log2FoldChange) > fc_in, padj < 0.05) %>%
      left_join(anno)
  #make result into excel csv
  write_csv(res_first, path = paste0("../Results/Kallisto_Gene/DE_GeneLevel_Kallisto_Results_Embryo_",fc_name,".csv"))

  #get result of heart mutant wt
  #add rowname for ensembl
  #filter FC and adjusted p val. abs mean no negative
  #add annotation
  res_second = as.data.frame(results(dds_fit, contrast = c("Var", "Heart_Mutant", "Heart_WT"))) %>%
       add_rownames("ensembl_gene_id") %>%
       filter(abs(log2FoldChange) > fc_in, padj < 0.05) %>%
       left_join(anno)
  #make result into excel csv
  write_csv(res_second, path = paste0("../Results/Kallisto_Gene//DE_GeneLevel_Kallisto_Results_Heart_",fc_name,".csv"))
}

#plot gene
#ensembl ID from anno
anno_in = anno %>% filter(zfin_id_symbol == "jupa")
#counts of jupa
#get jupa ensembl id from anno_in
#melt makes counts ok for ggplot
#attach pheno_in for colours
counts_df = as.data.frame(ncounts) %>% rownames_to_column("ensembl_gene_id") %>%
  filter(ensembl_gene_id == anno_in$ensembl_gene_id) %>% melt %>%
  left_join(pheno_in, by = c("variable" = "sample"))

#make plot
#group = 1 for colour to work
#facet like manual for looking at variables
#xlab ylab for labels
#title
graph = ggplot(counts_df, aes(variable,value, group=1)) +
               geom_point() + geom_line() +
               facet_grid(Tissue ~ Type, scales = "free_x") +
               theme(axis.text.x = element_blank()) +
               xlab("Samples") + ylab("Normalized Counts") +
               ggtitle("Gene Level JUPA")
print(graph)


library(pheatmap)
res_first = as.data.frame(results(dds_fit, contrast = c("Var", "Embryo_Mutant", "Embryo_WT"))) %>%
    add_rownames("ensembl_gene_id") %>%
    filter(abs(log2FoldChange) > fc_in, padj < 0.05) %>%
    left_join(anno)

#rlogTransformation for plot
#get differentially expressed genes
#get embryo sample, no heart
rld = assay(rlogTransformation(dds_fit))
rld = rld[res_first$ensembl_gene_id,]
rld = rld[,grepl("Embryo", pheno_in$Tissue)]
#make heatmap
pheatmap(rld,cluster_rows=TRUE,show_rownames=FALSE,cluster_cols=TRUE)


#goseq goterm kegg
#
library(goseq)

#get differentially expressed genes
#make DE col 1 is differential express 0 is not differential express
#use only DE col and ensembl colour
res_first = as.data.frame(results(dds_fit, contrast = c("Var", "Embryo_Mutant", "Embryo_WT"))) %>%
    add_rownames("ensembl_gene_id") %>%
    filter(abs(log2FoldChange) > fc_in, padj < 0.05) %>%
    mutate(DE = ifelse(abs(log2FoldChange) > log2(1.5) & padj < 0.05, 1,0),
           DE = ifelse(is.na(padj),0,DE))

#make goseq file, is 1 or 0
#name is ensembl gene id
goseq_de = res_first$DE
names(goseq_de) = res_first$ensembl_gene_id
#run from manual
pwf = nullp(goseq_de,"danRer6","ensGene", plot.fit = T)

#goterm enrichment from manual
#multiple test correction
#filter by p val
#write excel file
GO.wall = goseq(pwf,"danRer6","ensGene")
GO.wall$padj = p.adjust(GO.wall$over_represented_pvalue, method = "BH")
GO.wall.sig = GO.wall %>% filter(padj < 0.05) %>% arrange(padj)
write_csv(GO.wall.sig, path = "../Results/Kallisto_Gene/GOseq/GO_Terms_FC1.5.csv")

#kegg enrichment from manual
#multiple test correction
#filter by p val
#write excel file
GO.wall = goseq(pwf,"danRer6","ensGene",test.cats="KEGG")
GO.wall$padj = p.adjust(GO.wall$over_represented_pvalue,  method = "BH")
GO.wall.sig  = GO.wall %>% filter(padj < 0.05) %>% arrange(padj)
write_csv(GO.wall.sig, path = "../Results/Kallisto_Gene/GOseq/kegg_Terms_FC1.5.csv")
