# Simulation on effectiveness and robustness of multivariant cape model
# Parameters that user sets: n_v, n_p, n_i, allele_freq, int_list, ME_betas

library("MASS")
library("expm")
setwd("/Users/mingya/JAX/simulation/v1.3/bin")

source("helper_fn.R")


#############################################
### Define basic parameters:
#############################################
get.basic.params <- function() #fn only used for local scope
{
  source("tests/3v3p.R")
  
  return(list("n_v" = n_v, "n_p" = n_p, "n_i" = n_i, "allele_freq" = allele_freq, 
              "int_list" = int_list, "ME_betas" = ME_betas))
}

fabricate.data <- function(n_v, n_p, n_i, allele_freq, int_list, ME_betas)
{
  #############################################
  ### Create genetic population: X
  #############################################
  X <- create.genotype.population(n_v, allele_freq, n_i)
  genes <- X[,2:(n_v+1)]
  
  #############################################
  ### Define underlying genetic network: A
  #############################################
  colnames(int_list) <- c("v1", "v2", "weight")
  A <- get.var.adj.matrix(int_list, n_v)
  
  #############################################
  ### Calculate all betas using new model: betas
  #############################################
  betas <- get.betas.from.network(A, ME_betas)
  
  #############################################
  ### Calculate all phenotype values: Y
  #############################################
  Y <- X %*% betas
  
  #############################################
  ## Return objects
  #############################################
  return(list("X"=X, "Y"=Y, "A"=A, "betas"=betas))
}

#############################################
## Runtime: get params and fabricate data
#############################################
params <- get.basic.params()
dataset <- fabricate.data(params$n_v, params$n_p, params$n_i, params$allele_freq, 
                      params$int_list, params$ME_betas)
X <- dataset$X
Y <- dataset$Y
n_v <- params$n_v

save.image("fabricated_data.RData")
