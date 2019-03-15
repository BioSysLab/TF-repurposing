library(AnnotationDbi)
library(org.Hs.eg.db)
library(gdata)
library(dplyr)
library(tidyverse)
library(readr)
library(gplots)
library(xlsx)
library(RankAggreg)


###Analysis for Genetic Algorithm

#Load ranked lists and lists of scores
ordered_actives_up <- readRDS("ordered_myc_up.rds")
weights_actives_up <- readRDS("weights_myc_up.rds")

ordered_inactives <- readRDS("ordered_inactive.rds")
weights_inactives <- readRDS("weights_inactive.rds")


#Rank Aggregation with GA

GAW_act <- RankAggreg(ordered_actives_up, NCOL(ordered_actives_up) ,weights_actives_up, method= "GA",
                      distance="Spearman", seed=123, maxIter = 10000,
                      CP=.4, MP=.02, verbose=TRUE)

GAW_act_up <- GAW_act[["top.list"]]
final_act_up <- as.data.frame(GAW_act_up)
final_act_up <- final_act_up %>% mutate(rank=c(1:NCOL(ordered_actives_up)))
colnames(final_act_up) <- c("drugs","rank.up")
saveRDS(final_act_up,"final_up.rds")

#In downregulated we reverse the ranks in the end

GAW_act <- RankAggreg(ordered_inactives, NCOL(ordered_inactives) ,weights_inactives, method= "GA",
                      distance="Spearman", seed=123, maxIter = 10000,
                      CP=.4, MP=.02, verbose=TRUE)

GAW_act_down <- GAW_act[["top.list"]]
final_down <- as.data.frame(GAW_act_down)
final_down <- final_down %>% mutate(rank=c(1:NCOL(ordered_inactives)))
final_down %>% mutate(rank=max(rank)-rank+1)
final_down <- final_down[order(final_down$rank),]
colnames(final_down) <- c("drugs","rank.down")
saveRDS(final_down,"final_down.rds")


