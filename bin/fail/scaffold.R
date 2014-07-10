setwd("/Users/mingya/JAX/simulation/v1.4/bin/fail")
load("cape_little_cross.RData")
vv <- little.cross$var.to.var.p.val
vv = vv[vv[, "p.adjusted"] == 0, ]

#hist(little.cross$var.to.var.p.val[,"P_empirical"])
all.marks <- union(unique(vv[,1]), unique(vv[,2]))

clc <- function() cat(rep("\n", 50))
