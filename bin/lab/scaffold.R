setwd("~/JAX/simulation/v1.5/bin/lab/")

vv <- little.cross$var.to.var.p.val
vv = vv[vv[, "p.adjusted"] == 0, ]

#hist(little.cross$var.to.var.p.val[,"P_empirical"])
all.marks <- union(unique(vv[,1]), unique(vv[,2]))

clc <- function() cat(rep("\n", 50))