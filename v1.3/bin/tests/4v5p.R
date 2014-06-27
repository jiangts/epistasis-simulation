n_v <- 4 #num variant
n_p <- 5 #num pheno
n_i <- 1000 #num indiv
allele_freq <- 0.50 #50% allelic frequency.

int_list <- rbind(c(2,1,0.4), c(2,4,0.1), c(4,2,0.3), c(3,4,1), c(2,3,.07)) #list of interactions
ME_betas <- matrix(c(3, 3.5, 4, 3.5, 4, 3, 4, 3, 3.5, 1, 3, -2, -3, .4, 1, -3, -1, -2, 1, 2), 
                   nrow = n_v, ncol = n_p) #matrix of main effects