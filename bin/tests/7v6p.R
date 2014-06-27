library("reshape")

n_v <- 7 #num variant
n_p <- 6 #num pheno
n_i <- 1000 #num indiv
allele_freq <- 0.50 #50% allelic frequency.

int_list <- rbind(
  c(1,2,0.4),
  c(2,1,0.4),
  c(2,3,0.4),
  c(2,5,0.4),
  c(5,2,0.4),
  c(3,5,0.4),
  c(3,4,0.4),
  c(5,4,0.4),
  c(3,6,0.4),
  c(6,5,0.4),
  c(6,7,0.4)
  ) #list of interactions

set.seed(124)
ME_betas <- matrix(rnorm(n_v * n_p), 
                   nrow = n_v, ncol = n_p) #matrix of main effects