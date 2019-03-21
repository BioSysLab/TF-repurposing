library(fgsea)
library(gage)
library(EGSEAdata)
#library(KEGGREST)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
#library(cmapR)
#library(rhdf5)
#library(CARNIVAL)
#library(viper)
library(dplyr)
#data(kegg.gs)

#Load inputs----------------------------------------------------------------------------------
cmap_PA<- readRDS("cmap_PA.rds")
landmark <- read_tsv(file = "cmap_landmark_genes.txt")
egsea.data(species = "human",returnInfo = TRUE)
rownames(cmap_PA) <- landmark$`Entrez ID`
rownames(cmap_PA) <- as.character(rownames(cmap_PA))
cmap_myc_ranks_landmark<- readRDS("cmap_myc_ranks_landmark.rds")
cmap_myc_ranks_landmark <- cmap_myc_ranks_landmark %>% filter(myc_ranks<=18)
msigdb <- NULL
gsetdb.human <- NULL

##LOOP APPROACH-------------------------------------------------------------------------------
cmap_myc_ranks_landmark$pvalue.cellcycle <- 666
cmap_myc_ranks_landmark$pvalue.apopt <- 666
cmap_myc_ranks_landmark$padj.cellcycle <- 666
cmap_myc_ranks_landmark$padj.apopt <- 666
cmap_myc_ranks_landmark$rankES.cellcycle <- 666
cmap_myc_ranks_landmark$rankES.apopt <- 666
for (i in 1:NCOL(cmap_PA)) {
  paths <- fgsea(pathways = kegg.pathways$human$kg.sets, 
                    stats = cmap_PA[,i],
                    minSize=10,
                    maxSize=500,
                    nperm=10000)
  paths <- paths %>% mutate(ES.rank=rank(ES)) %>% select(pathway,pval,padj,ES.rank)
  cmap_myc_ranks_landmark$pvalue.cellcycle[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04110 Cell cycle") %>% select(pval)
  cmap_myc_ranks_landmark$pvalue.apopt[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04210 Apoptosis") %>% select(pval)
  cmap_myc_ranks_landmark$padj.apopt[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04210 Apoptosis") %>% select(padj)
  cmap_myc_ranks_landmark$padj.cellcycle[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04110 Cell cycle") %>% select(padj)
  cmap_myc_ranks_landmark$rankES.apopt[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04210 Apoptosis") %>% select(ES.rank)
  cmap_myc_ranks_landmark$rankES.cellcycle[which(cmap_myc_ranks_landmark$sig_id==colnames(cmap_PA)[i])]<-  paths %>% filter(pathway=="hsa04110 Cell cycle") %>% select(ES.rank)
  paths <- NULL
  print(paste("Iteration finished:",i))
  }

