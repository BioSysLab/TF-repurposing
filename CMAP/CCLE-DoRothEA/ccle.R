### Run DoRothEA on CCLE

library(cmapR)
library(rhdf5)
library(viper)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(grex)

### read ccle rna seq data
ds_path <- "CCLE_RNAseq_genes_rpkm_20180929.gct"

my_ds <- parse.gctx(ds_path)

ccle <- my_ds@mat

rownames(ccle) <- cleanid(rownames(ccle))

### Annotate the genes

ensembl <- rownames(ccle)
annoE <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = ensembl,
                               columns = c("SYMBOL", "GENENAME","ENTREZID"),
                               keytype = "ENSEMBL")

annoE <- annoE %>%
  filter(!is.na(SYMBOL))

annoE <- annoE[!duplicated(annoE$SYMBOL),]

ccle <- ccle[annoE$ENSEMBL,]

rownames(ccle) <- annoE$SYMBOL

### z-score transformation of expression

ccle_z <- t(scale(t(ccle)))

# Load TF regulon genesets in VIPER format
load("dorothea/TOP10score_viperRegulon.rdata")

# Clean TF names & explore object
names(viper_regulon) = sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)

# Run Dorothea-VIPER
TF_activities = viper(eset = ccle_z, regulon = viper_regulon, nes = T, method = 'none', minsize = 4, eset.filter = F)

### index of MYC
ind <- which(rownames(TF_activities) == "MYC_A")

### Ranks the TFs in the columns based on their score
TF_ranks <- apply(TF_activities,2,rank)

TF_ranks_myc <- TF_ranks[ind,]

saveRDS(TF_ranks_myc,"myc_ranked_cell_lines.rds")
