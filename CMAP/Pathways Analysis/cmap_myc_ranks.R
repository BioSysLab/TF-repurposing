library(tidyverse)
library(cmapR)
library(rhdf5)
library(CARNIVAL)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(viper)
library(readr)

files <- list.files(path = "cell_lines/",all.files = F,full.names = T,recursive = F)
ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"


### read landmark genes

landmark <- read_tsv(file = "C:/Users/user/Documents/deepCAGE/cmap/cmap_landmark_genes.txt")

### read dorothea module

load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL")) # loading the viper regulons

filler <- read_rds(files[1])
filler$myc_ranks <- 0
results <- filler[FALSE,]

for (i in 1:length(files)) {
  

  ### read the signature info for each cell line
  cell_sig <- read_rds(files[i])  
  
  ### get the cmap zscores for these signatures w/ landmark
  cmap <- get_cmap_signatures(cmap_path_to_gctx = ds_path,sig_ids = cell_sig$sig_id,landmark = TRUE,landmark_df = landmark)
  
  ### run viper
  viper <- runDoRothEA(cmap, regulon=viper_regulon, confidence_level=c('A','B','C'))
  
  ### rank the tfs for each drug
  viper_ranked <- apply(viper, MARGIN = 2, rank)
  
  ### index of myc
  ind <- which(rownames(viper_ranked) == "MYC")
  
  myc_ranks <- viper_ranked[ind,]
  
  myc_ranks <- as.data.frame(myc_ranks) %>%
    rownames_to_column("sig_id")
  
  cell_sig <- left_join(cell_sig,myc_ranks)
  
  results <- bind_rows(results,cell_sig)

  gc()
  print(i)
}

results <- results %>%
  group_by(pert_iname) %>%
  mutate(median_drug_ranks = median(myc_ranks)) %>%
  mutate(counts = n_distinct(cell_id)) %>%
  ungroup()


saveRDS(results,file = "cmap_myc_ranks_landmark.rds")

