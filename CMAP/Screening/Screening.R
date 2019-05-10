library(AnnotationDbi)
library(org.Hs.eg.db)
library(gdata)
library(dplyr)
library(tidyverse)
library(readr)
library(gplots)
library(xlsx)
library(GeneExpressionSignature)

#Drugs from Lopac Database
drugs_lopac <- read.xlsx("C:/Users/sofia/OneDrive/Documents/MyC/LOPAC_Compounds.xlsx",sheetIndex=1,as.data.frame=TRUE,
                         header=TRUE)
drugs_lopac <- drugs_lopac %>% select(Compound.Name)
colnames(drugs_lopac) <- c("drugs")
drugs_lopac$drugs <- toupper(drugs_lopac$drugs)
drugs_lopac <- drugs_lopac %>% mutate(exists_in_lopac="yes")

#Load drugs from our database
our_drugs <- readRDS("our_drugs.rds")

#Load 5 runs of pathway analysis
cmap.v1.1 <- readRDS("pathway_results_myc_cancerpath.rds")
cmap.v2.1 <- readRDS("pathway_results_myc_cancerpath2.0.rds")
cmap.v3.1 <- readRDS("pathway_results_myc_cancerpath3.0.rds")
cmap.v4.1 <- readRDS("pathway_results_myc_cancerpath4.0.rds")
cmap.v5.1 <- readRDS("pathway_results_myc_cancerpath5.0.rds")
cmap.v1.1$pert_iname <- toupper(cmap.v1.1$pert_iname) 
cmap.v2.1$pert_iname <- toupper(cmap.v2.1$pert_iname)
cmap.v3.1$pert_iname <- toupper(cmap.v3.1$pert_iname)
cmap.v4.1$pert_iname <- toupper(cmap.v4.1$pert_iname)
cmap.v5.1$pert_iname <- toupper(cmap.v5.1$pert_iname)

#First filter: Select drugs-cell lines pairs in which at least one of the desired pathways is activated
cmap.v1.1 <- cmap.v1.1 %>% filter(NES.cellcycle<0|NES.dnarep<0|NES.apopt>0|NES.natkiller>0|NES.pathcancer>0)
cmap.v2.1 <- cmap.v2.1 %>% filter(NES.cellcycle<0|NES.dnarep<0|NES.apopt>0|NES.natkiller>0|NES.pathcancer>0)
cmap.v3.1 <- cmap.v3.1 %>% filter(NES.cellcycle<0|NES.dnarep<0|NES.apopt>0|NES.natkiller>0|NES.pathcancer>0)
cmap.v4.1 <- cmap.v4.1 %>% filter(NES.cellcycle<0|NES.dnarep<0|NES.apopt>0|NES.natkiller>0|NES.pathcancer>0)
cmap.v5.1 <- cmap.v5.1 %>% filter(NES.cellcycle<0|NES.dnarep<0|NES.apopt>0|NES.natkiller>0|NES.pathcancer>0)

colnames(cmap.v1.1)[3] <- "drugs"
colnames(cmap.v2.1)[3] <- "drugs"
colnames(cmap.v3.1)[3] <- "drugs"
colnames(cmap.v4.1)[3] <- "drugs"
colnames(cmap.v5.1)[3] <- "drugs"

#Second Filter: Select the pairs in which at least one pathway is enriched and important
cmap.v1.2 <- cmap.v1.1 %>% filter((padj.cellcycle<0.05&rankNES.cellcycle<25)|(padj.apopt<0.05&rankNES.apopt>100)|(padj.dnarep<0.05&rankNES.dnarep<25)|(padj.natkiller<0.05&rankNES.natkiller>100)|(padj.pathcancer<0.05&rankNES.pathcancer>100))
cmap.v2.2 <- cmap.v2.1 %>% filter((padj.cellcycle<0.05&rankNES.cellcycle<25)|(padj.apopt<0.05&rankNES.apopt>100)|(padj.dnarep<0.05&rankNES.dnarep<25)|(padj.natkiller<0.05&rankNES.natkiller>100)|(padj.pathcancer<0.05&rankNES.pathcancer>100))
cmap.v3.2 <- cmap.v3.1 %>% filter((padj.cellcycle<0.05&rankNES.cellcycle<25)|(padj.apopt<0.05&rankNES.apopt>100)|(padj.dnarep<0.05&rankNES.dnarep<25)|(padj.natkiller<0.05&rankNES.natkiller>100)|(padj.pathcancer<0.05&rankNES.pathcancer>100))
cmap.v4.2 <- cmap.v4.1 %>% filter((padj.cellcycle<0.05&rankNES.cellcycle<25)|(padj.apopt<0.05&rankNES.apopt>100)|(padj.dnarep<0.05&rankNES.dnarep<25)|(padj.natkiller<0.05&rankNES.natkiller>100)|(padj.pathcancer<0.05&rankNES.pathcancer>100))
cmap.v5.2 <- cmap.v5.1 %>% filter((padj.cellcycle<0.05&rankNES.cellcycle<25)|(padj.apopt<0.05&rankNES.apopt>100)|(padj.dnarep<0.05&rankNES.dnarep<25)|(padj.natkiller<0.05&rankNES.natkiller>100)|(padj.pathcancer<0.05&rankNES.pathcancer>100))

#Get drugs derived from the data generated from each run of the pathway analysis
drugs_cmap1 <- cmap.v1.2 %>% select(drugs,median_drug_ranks)
drugs_cmap1 <- unique(drugs_cmap1)
colnames(drugs_cmap1) <- c("drugs","median_drug_ranks")
drugs_cmap1$median_drug_ranks <- drugs_cmap1$median_drug_ranks/171*100
drugs_cmap1 <- drugs_cmap1[order(drugs_cmap1$median_drug_ranks),]
drugs_cmap1 <- drugs_cmap1 %>% mutate(version=1)
saveRDS(drugs_cmap1,"drugs_cmap1.rds")
write.xlsx(as.data.frame(drugs_cmap1),"drugs_cmap1.xlsx", sheetName = "CMAP",col.names = TRUE, row.names = FALSE, append = FALSE)

drugs_cmap2 <- cmap.v2.2 %>% select(drugs,median_drug_ranks)
drugs_cmap2 <- unique(drugs_cmap2)
colnames(drugs_cmap2) <- c("drugs","median_drug_ranks")
drugs_cmap2$median_drug_ranks <- drugs_cmap2$median_drug_ranks/171*100
drugs_cmap2 <- drugs_cmap2[order(drugs_cmap2$median_drug_ranks),]
drugs_cmap2 <- drugs_cmap2 %>% mutate(version=2)
saveRDS(drugs_cmap2,"drugs_cmap2.rds")
write.xlsx(as.data.frame(drugs_cmap2),"drugs_cmap2.xlsx", sheetName = "CMAP",col.names = TRUE, row.names = FALSE, append = FALSE)

drugs_cmap3 <- cmap.v3.2 %>% select(drugs,median_drug_ranks)
drugs_cmap3 <- unique(drugs_cmap3)
colnames(drugs_cmap3) <- c("drugs","median_drug_ranks")
drugs_cmap3$median_drug_ranks <- drugs_cmap3$median_drug_ranks/171*100
drugs_cmap3 <- drugs_cmap3[order(drugs_cmap3$median_drug_ranks),]
drugs_cmap3 <- drugs_cmap3 %>% mutate(version=3)
saveRDS(drugs_cmap3,"drugs_cmap3.rds")
write.xlsx(as.data.frame(drugs_cmap3),"drugs_cmap3.xlsx", sheetName = "CMAP",col.names = TRUE, row.names = FALSE, append = FALSE)

drugs_cmap4 <- cmap.v4.2 %>% select(drugs,median_drug_ranks)
drugs_cmap4 <- unique(drugs_cmap4)
colnames(drugs_cmap4) <- c("drugs","median_drug_ranks")
drugs_cmap4$median_drug_ranks <- drugs_cmap4$median_drug_ranks/171*100
drugs_cmap4 <- drugs_cmap4[order(drugs_cmap1$median_drug_ranks),]
drugs_cmap4 <- drugs_cmap4 %>% mutate(version=4)
saveRDS(drugs_cmap4,"drugs_cmap4.rds")
write.xlsx(as.data.frame(drugs_cmap4),"drugs_cmap4.xlsx", sheetName = "CMAP",col.names = TRUE, row.names = FALSE, append = FALSE)

drugs_cmap5 <- cmap.v5.2 %>% select(drugs,median_drug_ranks)
drugs_cmap5 <- unique(drugs_cmap5)
colnames(drugs_cmap5) <- c("drugs","median_drug_ranks")
drugs_cmap5$median_drug_ranks <- drugs_cmap5$median_drug_ranks/171*100
drugs_cmap5 <- drugs_cmap5[order(drugs_cmap5$median_drug_ranks),]
drugs_cmap5 <- drugs_cmap5 %>% mutate(version=5)
saveRDS(drugs_cmap5,"drugs_cmap5.rds")
write.xlsx(as.data.frame(drugs_cmap5),"drugs_cmap5.xlsx", sheetName = "CMAP",col.names = TRUE, row.names = FALSE, append = FALSE)

#Combine and get only those drugs which appeared at least 3 times as final hits after the 2 filters
DRUGS <- rbind(drugs_cmap1,drugs_cmap2,drugs_cmap3,drugs_cmap4,drugs_cmap5)
DRUGS <- DRUGS %>% group_by(drugs) %>% summarize(version_counts=n_distinct(version))
final_drugs_cmpa <- DRUGS %>% filter(version_counts>=3)

cmap_top <- rbind(cmap.v1.2,cmap.v2.2,cmap.v3.2,cmap.v4.2,cmap.v5.2)
cmap_top <- left_join(cmap_top,DRUGS,by="drugs")
cmap_top <- cmap_top %>% filter(!is.na(version_counts)) %>% filter(version_counts>=3)
drugs_cmap_top <- cmap_top %>% select(drugs,median_drug_ranks)
drugs_cmap_top <- unique(drugs_cmap_top)
colnames(drugs_cmap_top) <- c("drugs","median_drug_ranks")
drugs_cmap_top$median_drug_ranks <- drugs_cmap_top$median_drug_ranks/171*100
drugs_cmap_top <- drugs_cmap_top[order(drugs_cmap_top$median_drug_ranks),]

#Save the suggested hits of the platform

saveRDS(drugs_cmap_top,"drugs_cmap_top.rds")
write.xlsx(as.data.frame(drugs_cmap_top),"drugs_cmap_top.xlsx", sheetName = "CMAP",col.names = TRUE, row.names = FALSE, append = FALSE)

##Find common drugs with lopac and ours, in order to quickly conduct experiments for validation

common2 <- left_join(cmap,drugs_lopac,by="drugs")
common2 <- common2 %>% filter(!is.na(exists_in_lopac))
saveRDS(common2,"common2_cancerpath.rds")
common3 <- left_join(cmap,our_drugs,by="drugs")
common3 <- common3 %>% filter(!is.na(ours))
saveRDS(common3,"common3_cancerpath.rds")


common2 <- left_join(cmap_top,drugs_lopac,by="drugs")
common2 <- common2 %>% filter(!is.na(exists_in_lopac))
common3 <- left_join(cmap_top,our_drugs,by="drugs")
common3 <- common3 %>% filter(!is.na(ours))
