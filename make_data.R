##################################################
# Making a dataset for analysis: genotype and phenotype information
# --> specify the following:
#   --> number of genes and phenotypes to create
#   --> genotypic distribution...
#   --> underlying genotype network
#   --> underlying main effects
#   --> model to use to get phenotype from genotype
#   --> noise level/noise type
########
# test performance of my two methodologies and cape(?)
##################################################
# for now, we  do 3 genes and 3 phenotypes
# there are 14 possible gene networks.
# we are doing this one:
#  A <--> B
#        /
#       v
#      C
##################################################
library("MASS") #get ginv() for generalized inverse
library("expm") #matrix exponentiation
setwd("/Users/mingya/JAX/SSP/bin/simulation/v1.1")
#set of helper fns that we use throughout
source("helper_fn.R")
#all the minimization code
source("direct_minimize.R")

#############################################
### Define basic parameters:
### number of variants, number of phenotypes,
### number of individuals in the dataset
#############################################
n_v <- 3 #num variant
n_p <- 3 #num pheno
n_i <- 1000 #num indiv

#############################################
### Define underlying genetic network: 
### set up the underlying variant interactions 
### that you want to test. Observe the format
### for interaction input. It will be converted
### into an adjacency matrix for you.
#############################################
AtoB <- c(1,2,0.9)
BtoA <- c(2,1,0.4)
BtoC <- c(2,3,0.8)
interactions <- rbind(AtoB, BtoA, BtoC)
colnames(interactions) <- c("v1", "v2", "weight")
int_adj <- get.var.adj.matrix(interactions, n_v) 
#^^ this is a delta matrix, NOT a m_ij matrix

#### the theoretical step!!! ####
### manipulating the adjacency matrix to create "overall effects"
#int_adj <- ginv(diag(n_v) - int_adj) - diag(n_v)
int_adj <- int_adj + int_adj %^% 2

#############################################
### Define the main effects:
### place these into a matrix. The columns
### are phenotypes, the rows are variants
### remember to put a leading row of 1's!
#############################################
# model to get phenotype values from genotype (and beta values)
main_effect_betas <- matrix(c(1, 3, 3.5, 4, 1, 3.5, 4, 3, 1, 4, 3, 3.5), 
                            nrow = n_p + 1, ncol = n_v)

### Calculating the full matrix of deltas,
### including main effect AND combined effect coefficients
betas <- get.beta.values(main_effect_betas, int_adj)

#############################################
### Create genotype distribution
### Assign a probability that a variant will
### exist at a locus.
#############################################
# create genotypic distribution
var_prob <- 0.50 #30% chance variant exists at locus.
design <- create.genotype.population(n_v, var_prob, n_i)

#############################################
### Create phenotype distribution
### Creates the dataset given the genotype 
### distribution, the underlying interactions,
### and the main effects.
#############################################
phenos <- design %*% betas

# TODO: incorporate noise of varying degrees...

########## "Solving" for our interaction values ##########
#cape method:
#############################################
### Using cape's methodology
### B_cape --> get coefficients of pairwise regression
###     matrix form
### delta_cape2 --> calculate deltas wrt 2 phenos
###     melted matrix form
### delta_cape --> calculate deltas wrt all phenos
###     melted matrix form
#############################################
B_cape <- ginv(design) %*% phenos #get betas
delta_cape2 <- get.cape.deltas2(B_cape, n_v)
delta_cape <- get.cape.deltas(B_cape, n_v)
print(signif(delta_cape, digits=3))

#separation method (get overall influnce, separate into direct and indirect)
#############################################
### Using the separation of overall influence method
### B_sep --> reshaped version of the delta_cape 
### melted matrix
### A_sep --> matrix of inferred direct influence
#############################################
B_sep <- get.sep.adj.matrix(delta_cape)
A_sep <- diag(n_v) - ginv(diag(n_v) + B_sep)
print(signif(A_sep, digits=3))


#direct search method (directly regress against many values of "delta")
#degree 2
#############################################
### Using the direct search of direct influence method
### result --> output of optim() function
### B_dir --> reshaping result into adjacency matrix
#############################################
result <- bfgs_opt()
B_dir <- melt(reshape.with.diag(result$par,n_v))

