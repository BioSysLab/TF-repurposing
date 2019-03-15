library(AnnotationDbi)
library(org.Hs.eg.db)
library(gdata)
library(dplyr)
library(tidyverse)
library(readr)
library(gplots)
library(xlsx)
library(RankAggreg)

#Final Top drugs

#Load merged ranked lists. In this case we load the lists from the GA.

#WARNING: For borda the final merged lists need reversing the ranks in inactives (in GA it is already done)

final_up<- readRDS("final_up.rds")
final_down<- readRDS("final_down.rds")


int <-  Reduce(intersect,list(final_up$drugs,
                                   final_down$drugs))
DrugScores_int <-  left_join(final_down,final_up, by="drugs")
DrugScores_int <-DrugScores_int %>% filter(!is.na(DrugScores_int$rank.down))
DrugScores_int <-DrugScores_int %>% filter(!is.na(DrugScores_int$rank.up))
DrugScores_int <- DrugScores_int %>% mutate(score=rank.up+rank.down)
DrugScores_int$score <- DrugScores_int$score/max(DrugScores_int$score)*100
DrugScores_int <- DrugScores_int[order(DrugScores_int$score),]
saveRDS(DrugScores_int,"drugs_score_borda.rds")
write.csv(DrugScores_int,"drugs_score_borda.csv")


