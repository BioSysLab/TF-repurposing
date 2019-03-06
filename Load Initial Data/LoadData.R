library(AnnotationDbi)
library(org.Hs.eg.db)
library(gdata)
library(dplyr)
library(tidyverse)
library(readr)
library(gplots)
library(xlsx)

#Load Datas
colu <- read.delim("v20._COLUMNS.txt")
ctrp_via <- read.delim("v20.data.curves_post_qc.txt")
sanger_via <- read.delim("v17.3_fitted_dose_response.txt",dec=",")
test <- as.numeric(sanger_via$LN_IC50)
sanger_via["LN_IC50"] <- exp(test)
toxdb <-NULL
compounds_ctrp <- read.delim("v20.meta.per_compound.txt") 
cell_lines_ctrp <- read.delim("v20.meta.per_cell_line.txt")
exp_id <- read.delim("v20.meta.per_experiment.txt")
exp1 <- left_join(exp_id,cell_lines_ctrp, by = "master_ccl_id")
exp2<- exp1[,c(1,2,4,7,9,10,13,14)]
ctrp <- left_join(compounds_ctrp,ctrp_via,by="master_cpd_id")
ctrp_final <- left_join(exp2,ctrp,by="experiment_id")
ctrp_final <- ctrp_final[,c(5,6,10,11,12,19,32)]


#Filter  and commons
sanger_filtered <- sanger_via
sanger_filtered <- sanger_filtered[,c(2,3,4,5,6,7,10)]
colnames(sanger_filtered)[3] <- "ccl_name"
colnames(sanger_filtered)[5] <- "cpd_name"
ctrp_filtered <-ctrp_final

sanger_filtered <- mutate(sanger_filtered,dataset="sanger")
sanger_filtered$cpd_name <- toupper(sanger_filtered$cpd_name)
ctrp_filtered <- mutate(ctrp_filtered,dataset="ctrp")
ctrp_filtered$cpd_name <-toupper(ctrp_filtered$cpd_name)


ctrp_duets <- unite(ctrp_filtered,"ccl/cpd",ccl_name,cpd_name,sep="/",remove=TRUE)
sanger_duets <- unite(sanger_filtered,"ccl/cpd",ccl_name,cpd_name,sep="/",remove=TRUE)
final_duets <- left_join(sanger_duets,ctrp_duets,by="ccl/cpd")
final_common_duets <- filter(final_duets,!is.na(final_duets$broad_cpd_id))
sanger_duets <- unite(sanger_filtered,"ccl/cpd",ccl_name,cpd_name,sep="/",remove=TRUE)
sanger_filtered2 <- filter(sanger_filtered,!(sanger_duets$`ccl/cpd` %in% final_common_duets$`ccl/cpd`))



#Combine

all_data <- dplyr::select(ctrp_filtered, 'ccl_name', 'cpd_name', IC50 = 'apparent_ec50_umol', 'dataset') %>% 
  rbind(., dplyr::select(sanger_filtered2, ccl_name, cpd_name, IC50 = 'LN_IC50', 'dataset'))
all_data <- all_data %>% mutate(dataset = as.factor(dataset))
all_data <- all_data %>% unite("ccl/cpd",ccl_name,cpd_name,sep="/",remove=FALSE)
all_data$ccl_name <- toupper(all_data$ccl_name)
all_data <- all_data %>% filter(!is.na(cpd_name))
all_data <- all_data %>% filter(!is.na(ccl_name))
all_data$IC50[which(all_data$IC50<0.000001)] <-0.000001
all_data$IC50[which(all_data$IC50>1000)] <-1000

saveRDS(all_data,file='all_data.rds')
