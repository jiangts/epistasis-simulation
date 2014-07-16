setwd("~/JAX/simulation/epistasis-simulation/bin/lab")
library("MASS")
load("steam_data.RData")

X<-little.cross$geno.for.pairscan * 2
Y<-little.cross$ET

#X.anal <- X[,c(25, 22, 5, 6, 7)]
#X.anal <- X[,c(18,23,76,83,75)]
#X.anal <- X[,c(3,37)]

#X.anal <- X[,c(32,33,50,51,52)]
X.anal <- X[,c(32,50,51)]
idxs = which(!is.na(rowSums(X.anal)))
X.anal <- X.anal[idxs,]
Y.anal <- Y[idxs,]

# need to define n_v
n_v <- ncol(X.anal)

# need to fix X to include x12, x123, etc.


for(k in 2:n_v)
{
  combos <- t(combn(c(1:n_v), k))
  combo_cols <- matrix(0, nrow=nrow(X.anal), ncol=nrow(combos))
  for(i in 1:nrow(combos))
  {
    combo_cols[,i] <- matrix(apply(X.anal[,combos[i,]], 1, prod)) #rowproducts
  }
  X.anal <- cbind(X.anal, combo_cols)
}

X.anal <- cbind(1, X.anal)


solve.betas <- ginv(X.anal) %*% Y.anal

source("helper_fn.R")

solve.bfgs <- bfgs()
solve.lev.marq <- lev.marq()
solve.nelder.mead <- nelder.mead()
solve.sim.anneal <- sim.anneal()

solve.bfgs.deltas <- reshape.with.diag(solve.bfgs$par, n_v)
solve.lev.marq.deltas <- reshape.with.diag(solve.lev.marq$par, n_v)
solve.nelder.mead.deltas <- reshape.with.diag(solve.nelder.mead$par, n_v)
solve.sim.anneal.deltas <- reshape.with.diag(solve.sim.anneal$par, n_v)

