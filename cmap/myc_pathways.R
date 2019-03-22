library(fgsea)
library(gage)
library(EGSEAdata)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)


#Load inputs----------------------------------------------------------------------------------
cmap_PA<- readRDS("C:/Users/sofia/OneDrive/Documents/cmap_PA.rds")
landmark <- read_tsv(file = "C:/Users/sofia/OneDrive/Documents/cmap_landmark_genes.txt")
egsea.data(species = "human",returnInfo = TRUE)
rownames(cmap_PA) <- landmark$`Entrez ID`
rownames(cmap_PA) <- as.character(rownames(cmap_PA))
cmap_myc_ranks_landmark<- readRDS("C:/Users/sofia/OneDrive/Documents/cmap_myc_ranks_landmark.rds")
cmap_myc_ranks_landmark <- cmap_myc_ranks_landmark %>% filter(myc_ranks<=18)
msigdb <- NULL
gsetdb.human <- NULL

##LOOP APPROACH-------------------------------------------------------------------------------
cmap_myc_ranks_landmark$pvalue.cellcycle <- 666
cmap_myc_ranks_landmark$pvalue.apopt <- 666
cmap_myc_ranks_landmark$pvalue.dnarep <- 666
cmap_myc_ranks_landmark$pvalue.natkiller <- 666
cmap_myc_ranks_landmark$padj.cellcycle <- 666
cmap_myc_ranks_landmark$padj.apopt <- 666
cmap_myc_ranks_landmark$padj.dnarep <- 666
cmap_myc_ranks_landmark$padj.natkiller <- 666
cmap_myc_ranks_landmark$rankNES.cellcycle <- 666
cmap_myc_ranks_landmark$rankNES.apopt <- 666
cmap_myc_ranks_landmark$rankNES.dnarep <- 666
cmap_myc_ranks_landmark$rankNES.natkiller <- 666
cmap_myc_ranks_landmark$NES.cellcycle <- 666
cmap_myc_ranks_landmark$NES.apopt <- 666
cmap_myc_ranks_landmark$NES.dnarep <- 666
cmap_myc_ranks_landmark$NES.natkiller <- 666
for (i in 1:NCOL(cmap_PA)) {
  paths <- fgsea(pathways = kegg.pathways$human$kg.sets, 
                 stats = cmap_PA[,i],
                 minSize=10,
                 maxSize=500,
                 nperm=10000)
  paths <- paths %>% mutate(NES.rank=rank(NES)) %>% select(pathway,pval,padj,NES.rank,NES)
  cmap_myc_ranks_landmark$pvalue.cellcycle[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04110 Cell cycle") %>% select(pval)
  cmap_myc_ranks_landmark$pvalue.apopt[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04210 Apoptosis") %>% select(pval)
  cmap_myc_ranks_landmark$pvalue.dnarep[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa03030 DNA replication") %>% select(pval)
  cmap_myc_ranks_landmark$pvalue.natkiller[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04650 Natural killer cell mediated cytotoxicity") %>% select(pval)
  
  cmap_myc_ranks_landmark$padj.apopt[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04210 Apoptosis") %>% select(padj)
  cmap_myc_ranks_landmark$padj.cellcycle[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04110 Cell cycle") %>% select(padj)
  cmap_myc_ranks_landmark$padj.dnarep[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa03030 DNA replication") %>% select(padj)
  cmap_myc_ranks_landmark$padj.natkiller[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04650 Natural killer cell mediated cytotoxicity") %>% select(padj)
  
  cmap_myc_ranks_landmark$rankNES.apopt[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04210 Apoptosis") %>% select(NES.rank)
  cmap_myc_ranks_landmark$rankNES.cellcycle[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04110 Cell cycle") %>% select(NES.rank)
  cmap_myc_ranks_landmark$rankNES.dnarep[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa03030 DNA replication") %>% select(NES.rank)
  cmap_myc_ranks_landmark$rankNES.natkiller[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04650 Natural killer cell mediated cytotoxicity") %>% select(NES.rank)
  
  cmap_myc_ranks_landmark$NES.apopt[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04210 Apoptosis") %>% select(NES)
  cmap_myc_ranks_landmark$NES.cellcycle[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04110 Cell cycle") %>% select(NES)
  cmap_myc_ranks_landmark$NES.dnarep[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa03030 DNA replication") %>% select(NES)
  cmap_myc_ranks_landmark$NES.natkiller[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04650 Natural killer cell mediated cytotoxicity") %>% select(NES)
  
  paths <- NULL
  print(paste("Iteration finished:",i))
}

