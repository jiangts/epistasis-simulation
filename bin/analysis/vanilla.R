setwd("/Users/mingya/JAX/simulation/v1.4/bin/analysis")

######################
# extract data
######################
#library("cape")
#lc = read.population("../../data/LitxLit_F2_cleaned.csv")
#X <- lc$geno
#Y <- lc$pheno
#detach("package:cape", unload=TRUE)
#save.image("little.cross.data.RData")

######################
# load data
######################
load("little.cross.data.RData")

filter.out.na <- function(X, Y){
  rs <- rowSums(is.na(X))
  #if more than 20 missing, mouse removed
  rs[rs>20] = NA #sum(!is.na(rs)) --> 1930 mice left
  X <- X[which(!is.na(rs)),]
  Y <- Y[which(!is.na(rs)),]
  
  cs <- colSums(is.na(X))
  #if more than 150 genotypes missing, genotype removed
  cs[cs>150] = NA #sum(!is.na(cs)) --> 75 genotypes left
  X <- X[,which(!is.na(cs))]
  
  Y <- Y[,5:22] #didn't take into account SEX
  rs <- rowSums(is.na(Y))
  rs[rs>5] = NA #sum(is.na(rs)) --> 1750 mice left
  
  return(list("X"=X, "Y"=Y))
}

filtered <- filter.out.na(X, Y)
X <- filtered$X
Y <- filtered$Y

