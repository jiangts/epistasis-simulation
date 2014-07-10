##Super trashy file!

setwd('/Users/mingya/JAX/simulation/v1.4/bin/comp')

# load data
library("cape")
rawdat <- read.population('../../data/LitxLit_F2_cleaned.csv')

dat <- delete.pheno(rawdat, c("masterID","mouseID2005", "pgm"))
dat <- create.covar(dat, "sex")
dat <- remove.ind.with.missing.pheno(dat) # unsure if needed...

#############################################################################
# concerns: how to choose the phenotypes? How to choose genotype "groups"?
# phenotype noise concerns. Number of individuals.
# Normalize phenotype data?
#############################################################################

#for(i in 1:ncol(dat$pheno))
#{
#  hist(dat$pheno[,i])
#} 
# everything is approximately normal, although Tb.N, BV.TV.TRABECULAR, 
# pctFat, body.weight don't look as "perfectly normal" as the others?


### CHOOSING PHENOTYPES Y ###
# black box-ish. Normalize, mean center, and SVD.
eigdat <- get.eigentraits(dat)
plotSVD(eigdat)
eigdat <- select.eigentraits(eigdat, traits.which = 1:5)


### CHOOSING GENOTYPES X ###
#1. to choose variants that will be included in the pair scan.
#2. large main eﬀects can obscure interactions. In suﬃciently powered studies, 
#conditioning on the large QTL can aid in the discovery of interactions, variants 
#with large main-eﬀects can be used as covariates in the pair scan.
eigdat <- singlescan(eigdat, n.perm=20, scan.what="eigentraits", 
                     auto.covar.selection = FALSE,verbose = TRUE)
save.image("eigdata.RData")

plotSinglescan(eigdat, mark.chr = TRUE, mark.covar = FALSE)

eigdat <- select.markers.for.pairscan(eigdat)
# generates geno.for.pairscan and covar.for.pairscan
#eigdat <- pairscan(eigdat, scan.what = "eigentraits", n.perm = 10, 
#                   min.per.genotype = 6, verbose = TRUE)
#save.image("eigdata-pair.RData")
#plotPairscan(eigdat, phenotype = c("ET1", "ET2", "ET3", "ET4", "ET5"), 
#             standardized = FALSE, pdf.label = "Pairscan.Regression.pdf")

eigdat <- error.prop(eigdat, perm = FALSE, verbose = TRUE)

X_all <- eigdat$geno.for.pairscan
X_all_mod <- X_all[,colSums(is.na(X_all)) < 400]
X_all_mod <- X_all_mod[rowSums(is.na(X_all_mod)) < 20,]
X<-X_all
Y <- eigdat$ET
X[is.na(X)]=0


##### STARTING FRESH HERE #####
little.cross$var.to.var.p.val[,"P_empirical"]
studentdata[studentdata$Drink == 'water',]

vv <- little.cross$var.to.var.p.val
vv1 <- subset(vv, vv[,"p.adjusted"] == 1)

block <- little.cross$linkage.blocks.collapsed
X <- little.cross$geno[,as.numeric(block)]
