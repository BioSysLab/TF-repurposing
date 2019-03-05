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

final_act_up<- readRDS("final_act_up.rds")
final_act_down<- readRDS("final_act_down.rds")
final_inact<- readRDS("final_inact.rds")

#Start merging and scoring

int_down <-  Reduce(intersect,list(final_inact$drugs,
                                   final_act_down$drugs))
DrugScores_int_down <-  left_join(final_inact,final_act_down, by="drugs")
DrugScores_int_down <-DrugScores_int_down %>% filter(!is.na(DrugScores_int_down$rank.inact))
DrugScores_int_down <-DrugScores_int_down %>% filter(!is.na(DrugScores_int_down$rank.act))
DrugScores_int_down <- DrugScores_int_down %>% mutate(score=rank.act+rank.inact)
saveRDS(DrugScores_int_down,"drugs_GAW_inact.rds")

int_up <-  Reduce(intersect,list(final_inact$drugs,
                                 final_act_down$drugs))
DrugScores_int_up <-  left_join(final_inact,final_act_up, by="drugs")
DrugScores_int_up <-DrugScores_int_up %>% filter(!is.na(DrugScores_int_up$rank.inact))
DrugScores_int_up <-DrugScores_int_up %>% filter(!is.na(DrugScores_int_up$rank.act))
DrugScores_int_up <- DrugScores_int_up %>% mutate(score=rank.act+rank.inact)
saveRDS(DrugScores_int_up,"drugs_GAW_inact.rds")

