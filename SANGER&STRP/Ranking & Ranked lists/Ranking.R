library(AnnotationDbi)
library(org.Hs.eg.db)
library(gdata)
library(dplyr)
library(tidyverse)
library(readr)
library(gplots)
library(xlsx)


#Sanger and CTRP DATA with IC50's. A quick data cleaning,organising etc
all_data <- readRDS("all_data.rds")
all_data_z <- all_data
all_data_z <- all_data_z %>% filter(!is.na(IC50))
all_data_z <- all_data_z %>% filter(!is.na(cpd_name))
all_data_z <- all_data_z %>% filter(!is.na(ccl_name))

#Only cell-lines with more than 10 different drugs tested on them remain on the data frame.
summary <- all_data_z %>% group_by(ccl_name) %>% summarise(n_distinct(cpd_name))
summary <- summary %>% filter(summary$`n_distinct(cpd_name)`>10)
all_data_z <- left_join(summary,all_data_z,by="ccl_name")
all_data_z <- all_data_z %>% filter(!is.na(`n_distinct(cpd_name)`))
all_data_z <- unique(all_data_z)

#Unique cell-lines in the data frame
cells <- as.data.frame(unique(all_data_z$ccl_name))
colnames(cells) <- "ccl_name"

#Count how many time a drugs is tested in a cell-line and separate those used once from those used more than once.
a <- all_data_z %>% group_by(ccl_name) %>% count(cpd_name)
a2 <- a %>% filter(n<2)
a3 <- a2 %>% unite("ccl/cpd",ccl_name,cpd_name,sep="/",remove=TRUE)
a4 <- a %>% filter(n>1)
a4 <- a4 %>% unite("ccl/cpd",ccl_name,cpd_name,sep="/",remove=TRUE)

#For those drugs tested more than once (probably with different dose,time etc) get the median of the IC50's from the experiment of the same drug in one cell-line
all_data_med <- all_data_z %>% filter(all_data_z$`ccl/cpd` %in% a4$`ccl/cpd`)
all_data_med <- all_data_med %>% group_by(`ccl/cpd`) %>% mutate(IC50=median(IC50))
all_data_med <- unique(all_data_med)

#For some pairs of cell-line/drug we have experiment both from sanger and ctrp and we choose those from ctrp
all_data_med <-all_data_med %>% group_by(`ccl/cpd`) %>% filter(dataset=='ctrp')

#Separate drugs with more IC50's values as explained above and merge the dataframes again in on one final data frame.
all_data_z <- as.data.frame(all_data_z %>% filter(all_data_z$`ccl/cpd` %in% a3$`ccl/cpd`))
all_data_test <- rbind.data.frame(all_data_z,all_data_med)
all_data_z <- all_data_test
all_data_test <- NULL
all_data_med <- NULL
a <- NULL
a2 <- NULL
a3 <- NULL
a4 <- NULL

#Load activity 

myc_ranked_cells <- readRDS("myc_ranked_cell_lines.rds")
myc_activity <- NULL
myc_activity$ccl_name <- NULL
myc_activity$mycrank <- NULL
myc_activity$ccl_name <- toupper(names(myc_ranked_cells))
myc_activity$mycrank <- myc_ranked_cells
myc_activity <- as.data.frame(myc_activity)
rownames(myc_activity) <- NULL
myc_activity <- myc_activity %>% separate(ccl_name, c("ccl_name", "tissue"), "_",remove = TRUE)

#threshold with window

myc_activity <- myc_activity %>% filter((mycrank<500 | mycrank>900))
myc_activity$activity <- 1
inactive_myc <- myc_activity %>% filter(mycrank<500)
inactive_myc$activity <- 0
active_myc <- myc_activity %>% filter(mycrank>900)
myc_activity <- rbind(active_myc,inactive_myc)


inactive_cellines <-myc_activity %>% filter(myc_activity$activity==0) %>% select(ccl_name,mycrank,activity)
active_cellines<-myc_activity %>% filter(myc_activity$activity==1) %>% select(ccl_name,mycrank,activity)
inactive_cellines <-arrange(inactive_cellines, desc(mycrank))
inactive_cellines <- inactive_cellines[(NROW(inactive_cellines)-NROW(active_cellines)+1):NROW(inactive_cellines),]
inactive_cellines <- unique(inactive_cellines)
myc_activity <- rbind(active_cellines,inactive_cellines)

###Ranking 

all_data_ranked <- all_data_z %>%group_by(ccl_name) %>%  mutate(rank=as.numeric(factor(rank(IC50))))
all_data_ranked <- left_join(all_data_ranked,myc_activity,by="ccl_name")
all_data_ranked <-  all_data_ranked %>% filter(!is.na(cpd_name))
all_data_ranked <-  all_data_ranked %>% filter(!is.na(ccl_name))
all_data_ranked <-  all_data_ranked %>% filter(!is.na(IC50))
all_data_ranked <-  all_data_ranked %>% filter(!is.na(activity))
actives <- all_data_ranked %>% filter(activity==1)
inactives <- filter(all_data_ranked,activity==0)
drugs_up <- as.data.frame(unique(actives$cpd_name))
colnames(drugs_up) <- "cpd_name"
drugs_down <- as.data.frame(unique(inactives$cpd_name))
colnames(drugs_down) <- "cpd_name"
drugs <- Reduce(intersect,list(drugs_down,drugs_up))

#Get ranked lists by also adding in every cell line the drugs that are missing


ordered_active_up <- actives %>% group_by(ccl_name) %>% arrange(rank,.by_group=TRUE) %>% select(cpd_name,`n_distinct(cpd_name)`,IC50,rank)
weights_active_up <- ordered_active_up %>% filter(`n_distinct(cpd_name)`>=600)
ordered_active_up <- ordered_active_up %>% filter(`n_distinct(cpd_name)`>=600)
ordered_active_up <- ordered_active_up %>% group_by(ccl_name) %>% mutate(medrank=median(rank),medIC50=median(IC50)) 
weights_active_up <- weights_active_up %>% group_by(ccl_name) %>% mutate(medrank=median(rank),medIC50=median(IC50)) 
weights_active_up <- ordered_active_up  %>% select(ccl_name,cpd_name,IC50,rank,medrank,medIC50)
ordered_active_up <- ordered_active_up  %>% select(ccl_name,cpd_name,IC50,rank,medrank,medIC50)
c1 <- unique(ordered_active_up$ccl_name)

#Add in every cell line the drugs that are missing by assuming they are in middle of the ranks
addition <- NULL
for (i in 1:length(c1)){
  tmp <- ordered_active_up %>% filter(ccl_name==c1[i])
  tmp <- left_join(drugs,tmp,by='cpd_name')
  tmp <- tmp %>% filter(is.na(IC50))
  tmp$ccl_name <- c1[i]
  tmp$`n_distinct(cpd_name)` <- NULL
  tmp$IC50 <- median(ordered_active_up$medIC50[which(ordered_active_up$ccl_name==c1[i])])
  tmp$rank <- median(ordered_active_up$medrank[which(ordered_active_up$ccl_name==c1[i])])
  tmp$medIC50 <- median(ordered_active_up$medIC50[which(ordered_active_up$ccl_name==c1[i])])
  tmp$medrank <- median(ordered_active_up$medrank[which(ordered_active_up$ccl_name==c1[i])])
  tmp <- tmp[c('ccl_name','cpd_name','IC50','rank','medrank','medIC50')]
  addition <- as.data.frame(rbind(addition,tmp))
}
ordered_active_up <- ordered_active_up %>% ungroup()
weights_active_up <-weights_active_up %>% ungroup()
ordered_active_up <- as.data.frame(rbind(ordered_active_up,addition))
weights_active_up <- as.data.frame(rbind(weights_active_up,addition))
ordered_active_up <- ordered_active_up %>% group_by(ccl_name) %>% arrange(rank,.by_group=TRUE) %>% select(ccl_name,cpd_name)
weights_active_up <- weights_active_up %>% group_by(ccl_name) %>% arrange(rank,.by_group=TRUE) %>% select(ccl_name,cpd_name,IC50)
addition <- NULL
tmp <- NULL

#Create matrix of weights and ranked lists with each row being a list
ordered_active_up1 <- NULL
ordered_active_up1 <- as.data.frame(rep(c(1:NROW(drugs)), times = length(c1)))
colnames(ordered_active_up1) <- "no"
ordered_active_up$no <-NULL
ordered_active_up$no <- ordered_active_up1$no
ordered_active_up <- ordered_active_up %>% spread(ccl_name,cpd_name)
ordered_active_up1 <- NULL
ordered_active_up$no <- NULL
ordered_active_up <- as.matrix(ordered_active_up)
ordered_active_up <- t(ordered_active_up)


weights_active_up <- weights_active_up %>% spread(cpd_name,IC50)
rownames(weights_active_up) <- weights_active_up$ccl_name
weights_active_up <- weights_active_up[,-1]
weights_active_up <- as.matrix(weights_active_up)

for (i in 1:NROW(weights_active_up)){
  weights_active_up[i,] <-weights_active_up[i,ordered_active_up[i,]] 
}
colnames(weights_active_up) <- NULL
rownames(weights_active_up) <- rownames(ordered_active_up)
weights_active_up <- 1/weights_active_up

saveRDS(ordered_active_up,"ordered_myc_up.rds")
saveRDS(weights_active_up,"weights_myc_up.rds")

#Myc Down ranked list. Same steps as used above in ranked active lists
ordered_inactives <- inactives %>% group_by(ccl_name) %>% arrange(rank,.by_group=TRUE) %>% select(cpd_name,`n_distinct(cpd_name)`,IC50,rank)
weights_inactives <- ordered_inactives %>% filter(`n_distinct(cpd_name)`>=600)
ordered_inactives <- ordered_inactives %>% filter(`n_distinct(cpd_name)`>=600) 
ordered_inactives <- ordered_inactives %>% group_by(ccl_name) %>% mutate(medrank=median(rank),medIC50=median(IC50)) 
weights_inactives <- weights_inactives %>% group_by(ccl_name) %>% mutate(medrank=median(rank),medIC50=median(IC50)) 
weights_inactives <- ordered_inactives  %>% select(ccl_name,cpd_name,IC50,rank,medrank,medIC50)
ordered_inactives <- ordered_inactives  %>% select(ccl_name,cpd_name,IC50,rank,medrank,medIC50)
c2 <- unique(ordered_inactives$ccl_name)

#Add in every cell line the drugs that are missing by assuming they are in middle of the ranks
addition <- NULL
for (i in 1:length(c2)){
  tmp <- ordered_inactives %>% filter(ccl_name==c2[i])
  tmp <- left_join(drugs,tmp,by='cpd_name')
  tmp <- tmp %>% filter(is.na(IC50))
  tmp$ccl_name <- c2[i]
  tmp$`n_distinct(cpd_name)` <- NULL
  tmp$IC50 <- median(ordered_inactives$medIC50[which(ordered_inactives$ccl_name==c2[i])])
  tmp$rank <- median(ordered_inactives$medrank[which(ordered_inactives$ccl_name==c2[i])])
  tmp$medIC50 <- median(ordered_inactives$medIC50[which(ordered_inactives$ccl_name==c2[i])])
  tmp$medrank <- median(ordered_inactives$medrank[which(ordered_inactives$ccl_name==c2[i])])
  tmp <- tmp[c('ccl_name','cpd_name','IC50','rank','medrank','medIC50')]
  addition <- as.data.frame(rbind(addition,tmp))
}
ordered_inactives <- ordered_inactives %>% ungroup()
weights_inactives <-weights_inactives %>% ungroup()
ordered_inactives <- as.data.frame(rbind(ordered_inactives,addition))
weights_inactives <- as.data.frame(rbind(weights_inactives,addition))
ordered_inactives <- ordered_inactives %>% group_by(ccl_name) %>% arrange(rank,.by_group=TRUE) %>% select(ccl_name,cpd_name)
weights_inactives <- weights_inactives %>% group_by(ccl_name) %>% arrange(rank,.by_group=TRUE) %>% select(ccl_name,cpd_name,IC50)
addition <- NULL
tmp <- NULL


ordered_inactives2 <- NULL
ordered_inactives2 <- as.data.frame(rep(c(1:NROW(drugs)), times = length(c2)))
colnames(ordered_inactives2) <- "no"
ordered_inactives$no <-NULL
ordered_inactives$no <- ordered_inactives2$no
ordered_inactives <- ordered_inactives %>% spread(ccl_name,cpd_name)
ordered_inactives2 <- NULL
ordered_inactives$no <- NULL
ordered_inactives <- as.matrix(ordered_inactives)
ordered_inactives <- t(ordered_inactives)

weights_inactives <- weights_inactives %>% spread(cpd_name,IC50)
rownames(weights_inactives) <- weights_inactives$ccl_name
weights_inactives <- weights_inactives[,-1]
weights_inactives <- as.matrix(weights_inactives)

for (i in 1:NROW(weights_inactives)){
  weights_inactives[i,] <-weights_inactives[i,ordered_inactives[i,]] 
}
colnames(weights_inactives) <- NULL
rownames(weights_inactives) <- rownames(ordered_inactives)
weights_inactives <- 1/weights_inactives

saveRDS(ordered_inactives,"ordered_inactive.rds")
saveRDS(weights_inactives,"weights_inactive.rds")

