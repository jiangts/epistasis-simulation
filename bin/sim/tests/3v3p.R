n_v <- 3 #num variant
n_p <- 3 #num pheno
n_i <- 1000 #num indiv
allele_freq <- 0.50 #50% allelic frequency.

int_list <- rbind(c(1,2,0.4), c(2,1,0.1), c(2,3,0.3)) #list of interactions
ME_betas <- matrix(c(3, 3.5, 4, 3.5, 4, 3, 4, 3, 3.5), 
                   nrow = n_v, ncol = n_p) #matrix of main effects