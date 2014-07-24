setwd("/Users/mingya/JAX/simulation/v1.5/bin/")
source("source.all.R")
source.all()

library(corpcor)
library(Matrix)
library(qpcR)
library(fdrtool)
library(evd)
library(igraph)
library(shape)

basedir <- "/home/bji/TCGA"
tag <- "subtype_ET4_0530_bycor"
#===============================================================
#specify the directories in which each class of data sits
#===============================================================
base.data.dir <- paste(basedir,"Data",sep="/") #a directory containing other data directories. This is where the data.obj will be saved after it is generated.
results.dir <- paste(basedir,"Results",tag,sep="/") #a directory where results will be stored

#############################

setwd(basedir)
if( is.na(match(paste("./", paste("Results",tag,sep="/") ,sep=""), list.dirs())) ) { dir.create( results.dir ) }

setwd(base.data.dir)
cross<- read.population("cape_cnv_05302014.csv")
summary(cross)

  
#================================================================================


  #covar=c("Basal","LumA", "Her2")
  covar=NULL
  covar.thresh = NULL
  pair.covar.pheno = c("pc1","pc2")

  use.pairs.thresh = FALSE
  max.pair.cor = 0.3 #only one of max.pair cor or min.per genotype can be used
  min.per.genotype = NULL
  scan.what = "rawtraits"
  transform.to.phenospace = FALSE
  pval.correction = "fdr"
  total.pair.perms = 1000000
  p.or.q = 0.001
  linkage.method = "genotype" #options are "genotype", "effects", or "prominence"
  #"genotype" requires a threshold for the minimum correlation
  #between linked markers.
  #"effects" groups makers by their main effects and interaction 
  #effects. It uses a hard threshold and a user-specified effect 
  #size drop to include only significant markers and adjacent markers 
  #within the effect size drop.
  #"prominence" finds all significant effects of each marker and
  #expands the linkage blocks based on the prominance of each peak
  #no additional parameters are used for this method.
  #all methods require a soft thresholding parameter for clustering
  #a cosine similarity matrix.
  #================================================================================
  
  setwd(results.dir)
  #cross <- create.covar(cross, covar)
  #cross <- norm.pheno(cross, mean.center = TRUE)
  #cross <- get.eigentraits(cross, scale.pheno = FALSE, normalize.pheno = FALSE)
  
  #pdf("svd.pdf")
  #plotSVD(cross, orientation = "vertical")
  #dev.off()
  
  #cross <- select.eigentraits(cross, traits.which = eig.which)
  
  cross<- singlescan(cross, n.perm = 3, covar = covar, scan.what=c("raw.traits"),alpha=c(0.01,0.05), verbose = TRUE)
  
  
  saveRDS(cross, "cross.RData")
  
  pdf("Singlescan_t_stat.pdf", width = 12)
  plotSinglescan(cross, mark.covar = TRUE, mark.chr = TRUE)

  dev.off()
  
  
  #cross <- set.covar(cross, pheno = pair.covar.pheno, markers = covar, is.covar = TRUE, plot.covar = FALSE)
  

  cross <- select.markers.for.pairscan(cross, use.pairs.threshold = FALSE)
  
  pdf("Selected.Markers.pdf", width = 12)
  plotSinglescan(cross, show.rejected.markers = TRUE)
  dev.off()
  
  #calculate the number of permutations to do per pair to 
  #achieve specified total number of permutations
  #num.pos.pairs <- choose(dim(cross$geno.for.pairscan)[2], 2)
  #pair.perms <- ceiling(total.pair.perms/num.pos.pairs)
  
  cross <- pairscan(cross, scan.what = scan.what, n.perm = 2, min.per.genotype = min.per.genotype, max.pair.cor = max.pair.cor, verbose = TRUE, num.pairs.limit = 10000000)
  
  saveRDS(cross, "cross.RData")
  
  plotPairscan(cross, phenotype = NULL, pdf.label = "Pairscan.Regression.pdf")
  
  
  cross <- error.prop(cross, perm = FALSE, verbose = TRUE)
  cross <- error.prop(cross, perm = TRUE, verbose = TRUE)
  
  plotPairscan(cross, phenotype = NULL, pdf.label = "test.Pairscan.Regression.pdf")
  saveRDS(cross, "cross.RData")
  
  cross <- calc.p(cross, pval.correction = pval.correction)
  plotPairscan(cross, phenotype = NULL, pdf.label = "test2.Pairscan.Regression.pdf")

  cross <- direct.influence(cross, transform.to.phenospace = transform.to.phenospace, verbose = TRUE, pval.correction = pval.correction, save.permutations = TRUE)
  
  
  saveRDS(cross, "cross.RData")
  
  pdf("subtype.variant.influences.pdf")
  plotVariantInfluences(cross, p.or.q = p.or.q, all.markers = FALSE, standardize = TRUE, not.tested.col = "lightgray", scale.effects=c("sqrt"),)
  dev.off()
  
  writeVariantInfluences(cross, p.or.q = max(c(p.or.q, 0.2)), filename = "Variant.Influences.csv")

  
  cross <- get.network(cross, p.or.q = p.or.q, collapse.linked.markers = FALSE)
  cross <- get.network(cross, p.or.q = p.or.q, collapse.linked.markers = TRUE, r.thresh=0.5)
  
  
  saveRDS(cross, "cross.RData")

  pdf("Network.Full.pdf", width = 12, height = 7)
  plotNetwork(cross, collapsed.net = FALSE)
  dev.off()
  
  pdf("Network.Collapsed.pdf", width = 12, height = 7)
  plotNetwork(cross, collapsed.net = TRUE)
  dev.off()	
  



