#### run borda for inactive

library(GeneExpressionSignature)

### read matrix of cell_lines x drugs with drug names inside

t11 <- readRDS("ordered_inactive.rds")

### create an embedding vector from drug names to integers

embedding <- as.data.frame(t(1:ncol(t11)))
colnames(embedding) <- t11[1,]
embedding <- as.matrix(embedding)

### reshape the matrix to drugs X cell_lines with ranks inside 

new_ranks <- matrix(data = 0,nrow = nrow(t11),ncol = ncol(t11))

### fill with ranks inside

for (i in 1:length(rownames(t11))) {
  
  new_ranks[i,] <- as.numeric(embedding[,t11[i,]])
  
}


rownames(new_ranks) <- rownames(t11)
colnames(new_ranks) <- colnames(embedding)
new_ranks <- t(new_ranks)

### create the dummy expression set to run borda

pheno <- as.matrix(rep("inactive",ncol(new_ranks)))
rownames(pheno) <- colnames(new_ranks)
pheno <- new("AnnotatedDataFrame",data=as.data.frame(pheno))
exprset <- new("ExpressionSet",exprs=new_ranks,phenoData=pheno)

### merge

merged <- RankMerging(exprset, MergingDistance = "Spearman", weighted = TRUE)

### result

merged_result <- exprs(merged)
### write results

write.csv(merged_result,"down.csv")
