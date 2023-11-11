
library(DESeq2)
library(tidyverse)
BiocManager::available()

setwd("/Users/bencekover/Library/CloudStorage/OneDrive-Personal/MSci Bahler lab/S.-Pombe-MLPs - Github")

map = read.csv("Bence folder/Analysis of Clement-Ziza RNA-seq/DESeq2/map.csv",row.names="sample")
map$snp = factor(map$snp)
raw_data = read.csv("Bence folder/Analysis of Clement-Ziza RNA-seq/DESeq2/raw_data.csv",row.names="X")



#safety check whether colnames in raw data are the same as in map
all(colnames(raw_data) == rownames(map))



#create deseq2 object

dds = DESeqDataSetFromMatrix(countData = raw_data,
                             colData = map,
                             design = ~ snp)


dds <- collapseReplicates(dds, dds$strain, dds$X)


#dds = dds[,keep]

dds$snp = relevel(dds$snp, ref = "0")

dds = DESeq(dds)

res = results(dds)


sig_res = results(dds, alpha=0.01)
sig_res= sig_res[order(sig_res$padj),]

summary(sig_res)

plotMA(sig_res)


install.packages("devtools")
devtools::install_github('kevinblighe/EnhancedVolcano')

library(EnhancedVolcano)


sig_res_final = sig_res[!grepl("SPNCRNA", rownames(sig_res)), ]


#short gene names
genes = read.delim("external data/Pombase files/gene_IDs_names_products.tsv", header=FALSE)

gene_names = list()
for (i in rownames(sig_res_final)) {
  if (i %in%genes$V1){
  gene_name = genes$V3[genes$V1 == i] 
  }
  else {
    gene_name = ""
  }

  if (gene_name == ""){
    gene_name = "" }

  gene_names = append(gene_names,gene_name) 
  }
  




EnhancedVolcano(sig_res_final,
                lab = gene_names,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01/length(raw_data),
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4.0,
                drawConnectors=TRUE)


write.csv(sig_res_final, file="Bence folder/Analysis of Clement-Ziza RNA-seq/DESeq2/DE_results.csv")


#The final plot for the publication was actually made in Python.


