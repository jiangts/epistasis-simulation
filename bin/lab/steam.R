setwd("~/JAX/simulation/v1.5/bin/lab/")
library("MASS")
source("helper_fn.R")
load("steam_data.RData")

X<-little.cross$geno.for.pairscan #* 2
Y<-little.cross$ET

data = get.relevant.markers(c(32,50,51))
X = data$X
Y = data$Y

# need to define n_v
n_v <- ncol(X)

design = get.design.mat(X,n_v)

solve.betas <- ginv(design) %*% Y
solve.all.betas = solve.betas

solve.lev.marq <- lev.marq()

CV.dist = generate.CV.distribution(design, maxerr=.002, cv.sample.size = 101)

p.values = get.p.values(CV.dist)
