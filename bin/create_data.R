# Simulation on effectiveness and robustness of multivariant cape model
# Parameters that user sets: n_v, n_p, n_i, allele_freq, int_list, ME_betas

library("MASS")
library("expm")
setwd("/Users/mingya/JAX/simulation/v1.3/bin")

source("helper_fn.R")

#############################################
### Define basic parameters:
#############################################
n_v <- 3 #num variant
n_p <- 3 #num pheno
n_i <- 1000 #num indiv

#############################################
### Create genetic population: X
#############################################
allele_freq <- 0.50 #50% allelic frequency.
X <- create.genotype.population(n_v, allele_freq, n_i)
genes <- X[,2:(n_v+1)]

#############################################
### Define underlying genetic network: A
#############################################
int_list <- rbind(c(1,2,0.4), c(2,1,0.1), c(2,3,0.3))
colnames(int_list) <- c("v1", "v2", "weight")
A <- get.var.adj.matrix(int_list, n_v)

#############################################
### Define main effects: ME
#############################################
ME_betas <- matrix(c(3, 3.5, 4, 3.5, 4, 3, 4, 3, 3.5), 
                   nrow = n_p, ncol = n_v)

#############################################
### Calculate all betas using new model: betas
#############################################
betas <- get.betas.from.network(A, ME_betas)

#############################################
### Calculate all phenotype values: Y
#############################################
Y <- X %*% betas

#############################################
## Save workspace
#############################################
save.image("fabricated_data.RData")