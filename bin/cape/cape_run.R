setwd("~/JAX/simulation/v1.5/bin/cape")

#==================================================
# load the cape library and the example data
#==================================================
library("cape")
little.cross <- read.population('../../data/LitxLit_F2_cleaned.csv')


#==================================================
# look at the structure of the data
#==================================================
str(little.cross)



#==================================================
# select the desired phenotypes for analysis
#==================================================
little.cross <- delete.pheno(little.cross, 
                             phenotypes = c("masterID","mouseID2005", "pgm", "Final.IGF.1"))



#==================================================
# specify the variable "mom" as a covariate.
# this function moves the variable from the 
# phenotype matrix to the genotype matrix
#==================================================
little.cross <- create.covar(little.cross, c("sex"))



#==================================================
# plot the distributions of the phenotypes. There
# is no single function to do this in CAPE yet, 
# but one is coming.
#==================================================
#layout(matrix(c(1:3), nrow = 1))
#hist.data <- apply(matrix(c(1:dim(little.cross$pheno)[2]), nrow = 1), 2, function(x) 
#  hist(little.cross$pheno[,x], main = colnames(little.cross$pheno)[x], 
#       xlab = colnames(little.cross$pheno)[x]))


#==================================================
# mean-center the phenotypes and transform them to
# fit a normal distribution
#==================================================
little.cross <- norm.pheno(little.cross, mean.center = TRUE)

#==================================================
# plot the distributions of the phenotypes again to
# make sure they were normalized
#==================================================
#layout(matrix(c(1:3), nrow = 1))
#hist.data <- apply(matrix(c(1:dim(little.cross$pheno)[2]), nrow = 1), 2, function(x) 
#  hist(little.cross$pheno[,x], main = colnames(little.cross$pheno)[x], 
#       xlab = colnames(little.cross$pheno)[x]))



#==================================================
# plot the correlations of the normalized phenotypes
# no CAPE function exists for doing this yet, but
# one is on the way
#==================================================
#layout(matrix(c(1:3), nrow = 1))
#i <- 1
#while(i < dim(little.cross$pheno)[2]){
#  j <- i + 1
#  while(j <= dim(little.cross$pheno)[2]){
#    plot(little.cross$pheno[,i], little.cross$pheno[,j], 
#         xlab = colnames(little.cross$pheno)[i], ylab = 
#           colnames(little.cross$pheno)[j], main = paste("r =", 
#                                                          round(cor(little.cross$pheno[,i], little.cross$pheno[,j], 
#                                                                    use = "complete.obs"), 2)), cex.lab = 1.2)
#    j <- j + 1
#  }
#  i = i + 1
#}



#==================================================
# use singular value decomposition to calculate 
# eigentraits (ET).
#==================================================
little.cross <- get.eigentraits(little.cross, scale.pheno = FALSE,
                                 normalize.pheno = FALSE)



#==================================================
# plot the result of the SVD. This will help in
# selecting which ET will be used for the
# analysis.
#==================================================
#plotSVD(little.cross, orientation = "vertical")



#==================================================
# select which ET will be used in the analysis
# here we select the first and second
#==================================================
little.cross <- select.eigentraits(little.cross, traits.which = 1:5)



#==================================================
# perform linear regression to associate each marker
# with each ET. We won't use alpha.for.covar or
# alpha.for.pairs to select markers, but these values
# will show up on the plot, so it's nice to set them
# at values you care about.
#==================================================
little.cross <- singlescan(little.cross, n.perm = 10, covar = c("sex"), 
                            scan.what = "eigentraits", alpha = c(0.01, 0.05),
                            verbose = TRUE)



#==================================================
# plot the results of the regression
#==================================================
plotSinglescan(little.cross)



#==================================================
# make sure "mom" is classified as a covariate for
# the pairwise marker scan. If there are other 
# markers you want to include as covariates, they
# can be specified here.
#==================================================
little.cross <- set.covar(little.cross, pheno = "ET1",  
                           markers = c("sex"), is.covar = TRUE, plot.covar = TRUE)


#==================================================
# Select markers for the pairscan. This step must
# be done even if no threhold is used. It checks
# for identical pairs of markers in addition to
# filtering.
#==================================================
little.cross <- select.markers.for.pairscan(little.cross, 
                                             use.pairs.threshold = FALSE)


#==================================================
# Perform the pairwise regression. This associates
# each marker pair with each ET and collects all 
# the coefficients from each regression model.
#==================================================
little.cross <- pairscan(little.cross, scan.what = "eigentraits", 
                          n.perm = 3, min.per.genotype = 6, verbose = TRUE)



#==================================================
# plot the interaction coefficients from the 
# pairwise regression
#==================================================
plotPairscan(little.cross, standardized = TRUE)



#==================================================
# Perform the error propagation step. This ensures
# that huge effects aren't counted as significant
# if they also have huge standard errors.
#==================================================
little.cross <- error.prop(little.cross, perm = FALSE, verbose = TRUE)


#==================================================
# perform the error propagation step on all the
# permuted data as well. This is the longest 
# step of the whole procedure.
#==================================================
little.cross <- error.prop(little.cross, perm = TRUE, verbose = TRUE)
#UGH


#==================================================
# calculate a p value for each genetic interaction
# choose the method of correction for multiple testing
# "fdr" is less stringent than "holm." "holm" is 
# essentially a bonferoni correction.
#==================================================
little.cross <- calc.p(little.cross, pval.correction = "holm")


#==================================================
# calculate the main effect of each marker on each
# phenotype.
#==================================================
little.cross <- direct.influence(little.cross, transform.to.phenospace = TRUE, pval.correction = "holm")



#==================================================
# plot the results as an asymmetric matrix showing
# source and target markers
#==================================================
plotVariantInfluences(little.cross, p.or.q = 0.05, 
                      all.markers = FALSE, standardize = FALSE, not.tested.col = "lightgray")


#==================================================
# write a table of the of the variant-to-variant
# influences
#==================================================
writeVariantInfluences(little.cross, p.or.q = 0.05, 
                       filename =  "Variant.Influences.csv", delim = ",",
                       mark.covar = FALSE)



#==================================================
# calculate the network from the adjacency matrix
# calculate both the full network and a network
# in which the loci are collapsed into QTL based on
# linkage between markers.
#==================================================
little.cross <- get.network(little.cross, p.or.q = 0.05, collapse.linked.markers = FALSE)
little.cross <- get.network(little.cross, p.or.q = 0.05, collapse.linked.markers = TRUE)


#==================================================
# plot the networks
#==================================================
plotNetwork(little.cross, collapsed.net = FALSE)

#==================================================
# Unfortunately, there seems to be a bug in plotting
# the collapsed network right now. I will fix this
# as soon as possible and put up a working version
# to use with your own data.
#==================================================
plotNetwork(little.cross, collapsed.net = TRUE)

