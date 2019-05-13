library(tidyverse)
library(rcdk)
library(rJava)
library(iterators)
#PARSING FOR LOPAC
iter <- iload.molecules(file.path("data/LOPAC/LOPAC.sdf"), type='sdf',skip = TRUE)
lopac_mol <- list()

i <- 0
while(hasNext(iter)) {
  i = i+1
  lopac_mol[[i]] <- nextElem(iter)
  
}

lopac_names <- lapply(lopac_mol, get.properties)
lopac_names_only <- lopac_names %>% lapply(. %>% .$ID) %>% unlist
names(lopac_mol) <- lopac_names_only


### LOPAC smiles

lopac_smiles <- lapply(lopac_mol,get.smiles, flavor = smiles.flavors(c('Generic')))

lopac_smiles2 <- as.data.frame(unlist(lopac_smiles))

lopac_smiles2 <- lopac_smiles2 %>%
  rownames_to_column("ID")

colnames(lopac_smiles2) <- c("ID","smiles")

### save the smiles of lopac

saveRDS(lopac_smiles2,"data/LOPAC/lopac_smiles.rds")
