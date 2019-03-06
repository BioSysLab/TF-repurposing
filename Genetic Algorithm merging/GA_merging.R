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

#Load active ranked lists and lists of scores
ordered_actives_up <- readRDS("ordered_active_up.rds")
weights_actives_up <- readRDS("weights_active_up.rds")

ordered_actives_down <- readRDS("ordered_active_down.rds")
weights_actives_down <- readRDS("weights_active_down.rds")

ordered_inactives <- readRDS("ordered_inactive_downanalysis.rds")
weights_inactives <- readRDS("weights_inactive_downanalysis.rds")


#Rank Aggregation with GA

GAW_act <- RankAggreg(ordered_actives_up, NCOL(ordered_actives_up) ,weights_actives_up, method= "GA",
                      distance="Spearman", seed=123, maxIter = 10000,
                      CP=.4, MP=.02, verbose=TRUE)

GAW_act_up <- GAW_act[["top.list"]]
final_act_up <- as.data.frame(GAW_act_up)
final_act_up <- final_act_up %>% mutate(rank=c(1:NCOL(ordered_actives_up)))
colnames(final_act_up) <- c("drugs","rank.act")
saveRDS(final_act_up,"final_act_up.rds")

GAW_act <- RankAggreg(ordered_actives_down, NCOL(ordered_actives_down) ,weights_actives_down, method= "GA",
                      distance="Spearman", seed=123, maxIter = 10000,
                      CP=.4, MP=.02, verbose=TRUE)

GAW_act_down <- GAW_act[["top.list"]]
final_act_down <- as.data.frame(GAW_act_down)
final_act_down <- final_act_down %>% mutate(rank=c(1:NCOL(ordered_actives_down)))
colnames(final_act_down) <- c("drugs","rank.act")
saveRDS(final_act_down,"final_act_down.rds")

#In inactives we reverse the ranks in the end

GAW_inact <- RankAggreg(ordered_inactives, NCOL(ordered_inactives), weights_inactives, method= "GA",
                        distance="Spearman", seed=123, maxIter = 10000,
                        CP=.4, MP=.02, verbose=TRUE)

GAW_inact <- GAW_inact[["top.list"]]
final_inact <- as.data.frame(GAW_inact)
final_inact <- final_inact %>% mutate(rank=c(1:NCOL(ordered_inactives)))
final_inact <- final_inact %>% mutate(rank=max(rank)-rank+1)
final_inact <- final_inact[order(final_inact$rank),]
colnames(final_inact) <- c("drugs","rank.inact")
saveRDS(final_inact,"final_inact.rds")

