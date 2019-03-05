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

ordered_iactives <- readRDS("ordered_inactive_downanalysis.rds")
weights_inactives <- readRDS("weights_inactive_downanalysis.rds")


#Rank Aggregation with GA

GAW_act <- RankAggreg(ordered_actives_up, 319 ,weights_actives_up, method= "GA",
                      distance="Spearman", seed=123, maxIter = 10000,
                      CP=.4, MP=.02, verbose=TRUE)

GAW_act_up <- GAW_act[["top.list"]]
final_act_up <- as.data.frame(GAW_act_up)
final_act_up <- final_act_up %>% mutate(rank=c(1:319))
colnames(final_act_up) <- c("drugs","rank.act")
saveRDS(final_act_up,"final_act_up.rds")

GAW_act <- RankAggreg(ordered_actives_down, 448 ,weights_actives_down, method= "GA",
                      distance="Spearman", seed=123, maxIter = 10000,
                      CP=.4, MP=.02, verbose=TRUE)

GAW_act_down <- GAW_act[["top.list"]]
final_act_down <- as.data.frame(GAW_act_down)
final_act_down <- final_act_down %>% mutate(rank=c(1:448))
colnames(final_act_down) <- c("drugs","rank.act")
saveRDS(final_act_down,"final_act_down.rds")

GAW_inact <- RankAggreg(ordered_inactives, 435, weights_inactives, method= "GA",
                        distance="Spearman", seed=123, maxIter = 10000,
                        CP=.4, MP=.02, verbose=TRUE)

GAW_inact <- GAW_inact[["top.list"]]
final_inact <- as.data.frame(GAW_inact)
final_inact <- final_inact %>% mutate(rank=c(1:435))
final_inact <- final_inact %>% mutate(rank=max(rank)-rank+1)
final_inact <- final_inact[order(final_inact$rank),]
colnames(final_inact) <- c("drugs","rank.inact")
saveRDS(final_inact,"final_inact.rds")

