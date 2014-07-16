setwd("~/JAX/simulation/epistasis-simulation/bin/lab")

#==================================================
# load the cape library and the example data
#==================================================
library("cape")
obesity.cross <- read.population('../../data/LitxLit_F2_cleaned.csv')


#==================================================
# look at the structure of the data
#==================================================
str(obesity.cross)



#==================================================
# select the desired phenotypes for analysis
#==================================================
obesity.cross <- delete.pheno(obesity.cross, 
                              phenotypes = c("masterID","mouseID2005", "pgm"))



#==================================================
# specify the variable "mom" as a covariate.
# this function moves the variable from the 
# phenotype matrix to the genotype matrix
#==================================================
obesity.cross <- create.covar(obesity.cross, c("sex", "Final.IGF.1"))



#==================================================
# plot the distributions of the phenotypes. There
# is no single function to do this in CAPE yet, 
# but one is coming.
#==================================================
#layout(matrix(c(1:3), nrow = 1))
#hist.data <- apply(matrix(c(1:dim(obesity.cross$pheno)[2]), nrow = 1), 2, function(x) 
#  hist(obesity.cross$pheno[,x], main = colnames(obesity.cross$pheno)[x], 
#       xlab = colnames(obesity.cross$pheno)[x]))


#==================================================
# mean-center the phenotypes and transform them to
# fit a normal distribution
#==================================================
obesity.cross <- norm.pheno(obesity.cross, mean.center = TRUE)

#==================================================
# plot the distributions of the phenotypes again to
# make sure they were normalized
#==================================================
#layout(matrix(c(1:3), nrow = 1))
#hist.data <- apply(matrix(c(1:dim(obesity.cross$pheno)[2]), nrow = 1), 2, function(x) 
#  hist(obesity.cross$pheno[,x], main = colnames(obesity.cross$pheno)[x], 
#       xlab = colnames(obesity.cross$pheno)[x]))



#==================================================
# plot the correlations of the normalized phenotypes
# no CAPE function exists for doing this yet, but
# one is on the way
#==================================================
#layout(matrix(c(1:3), nrow = 1))
#i <- 1
#while(i < dim(obesity.cross$pheno)[2]){
#  j <- i + 1
#  while(j <= dim(obesity.cross$pheno)[2]){
#    plot(obesity.cross$pheno[,i], obesity.cross$pheno[,j], 
#         xlab = colnames(obesity.cross$pheno)[i], ylab = 
#           colnames(obesity.cross$pheno)[j], main = paste("r =", 
#                                                          round(cor(obesity.cross$pheno[,i], obesity.cross$pheno[,j], 
#                                                                    use = "complete.obs"), 2)), cex.lab = 1.2)
#    j <- j + 1
#  }
#  i = i + 1
#}



#==================================================
# use singular value decomposition to calculate 
# eigentraits (ET).
#==================================================
obesity.cross <- get.eigentraits(obesity.cross, scale.pheno = FALSE,
                                 normalize.pheno = FALSE)



#==================================================
# plot the result of the SVD. This will help in
# selecting which ET will be used for the
# analysis.
#==================================================
#plotSVD(obesity.cross, orientation = "vertical")



#==================================================
# select which ET will be used in the analysis
# here we select the first and second
#==================================================
obesity.cross <- select.eigentraits(obesity.cross, traits.which = 1:5)



#==================================================
# perform linear regression to associate each marker
# with each ET. We won't use alpha.for.covar or
# alpha.for.pairs to select markers, but these values
# will show up on the plot, so it's nice to set them
# at values you care about.
#==================================================
obesity.cross <- singlescan(obesity.cross, n.perm = 10, covar = c("sex","Final.IGF.1"), 
                            scan.what = "eigentraits", auto.covar.selection = FALSE, alpha.for.covar = 0.01, 
                            alpha.for.pairs = 0.05, verbose = TRUE)



#==================================================
# plot the results of the regression
#==================================================
plotSinglescan(obesity.cross)



#==================================================
# make sure "mom" is classified as a covariate for
# the pairwise marker scan. If there are other 
# markers you want to include as covariates, they
# can be specified here.
#==================================================
obesity.cross <- set.covar(obesity.cross, pheno = "ET1",  
                           markers = c("sex", "Final.IGF.1"), is.covar = TRUE, plot.covar = TRUE)


#==================================================
# Select markers for the pairscan. This step must
# be done even if no threhold is used. It checks
# for identical pairs of markers in addition to
# filtering.
#==================================================
obesity.cross <- select.markers.for.pairscan(obesity.cross, 
                                             use.pairs.threshold = FALSE)


#==================================================
# Perform the pairwise regression. This associates
# each marker pair with each ET and collects all 
# the coefficients from each regression model.
#==================================================
obesity.cross <- pairscan(obesity.cross, scan.what = "eigentraits", 
                          n.perm = 3, min.per.genotype = 6, verbose = TRUE)



#==================================================
# plot the interaction coefficients from the 
# pairwise regression
#==================================================
plotPairscan(obesity.cross, standardized = TRUE)



#==================================================
# Perform the error propagation step. This ensures
# that huge effects aren't counted as significant
# if they also have huge standard errors.
#==================================================
obesity.cross <- error.prop(obesity.cross, perm = FALSE, verbose = TRUE)


#==================================================
# perform the error propagation step on all the
# permuted data as well. This is the longest 
# step of the whole procedure.
#==================================================
obesity.cross <- error.prop(obesity.cross, perm = TRUE, verbose = TRUE)
#UGH


#==================================================
# calculate a p value for each genetic interaction
# choose the method of correction for multiple testing
# "fdr" is less stringent than "holm." "holm" is 
# essentially a bonferoni correction.
#==================================================
obesity.cross <- calc.p(obesity.cross, pval.correction = "holm")


#==================================================
# calculate the main effect of each marker on each
# phenotype.
#==================================================
obesity.cross <- direct.influence(obesity.cross, transform.to.phenospace = TRUE, pval.correction = "holm")



#==================================================
# plot the results as an asymmetric matrix showing
# source and target markers
#==================================================
plotVariantInfluences(obesity.cross, p.or.q = 0.05, 
                      all.markers = FALSE, standardize = FALSE, not.tested.col = "lightgray")


#==================================================
# write a table of the of the variant-to-variant
# influences
#==================================================
writeVariantInfluences(obesity.cross, p.or.q = 0.05, 
                       filename =  "Variant.Influences.csv", delim = ",",
                       mark.covar = FALSE)



#==================================================
# calculate the network from the adjacency matrix
# calculate both the full network and a network
# in which the loci are collapsed into QTL based on
# linkage between markers.
#==================================================
obesity.cross <- get.network(obesity.cross, p.or.q = 0.05, collapse.linked.markers = FALSE)
obesity.cross <- get.network(obesity.cross, p.or.q = 0.05, collapse.linked.markers = TRUE)


#==================================================
# plot the networks
#==================================================
plotNetwork(obesity.cross, collapsed.net = FALSE)

#==================================================
# Unfortunately, there seems to be a bug in plotting
# the collapsed network right now. I will fix this
# as soon as possible and put up a working version
# to use with your own data.
#==================================================
plotNetwork(obesity.cross, collapsed.net = TRUE)

