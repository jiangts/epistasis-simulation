setwd("/Users/mingya/JAX/simulation/v1.4/bin")

source("create_data.R")

#############################################
### Perform linear regression to get beta matrix
#############################################
solve.betas <- ginv(X) %*% Y
#solve.deltas <- nonlinear.squares.regression

# note on code optimization: The slow stuff is everything touched by cost_fn,
# which will be run many many times. Look there for speed improvements
solve.bfgs <- bfgs()
solve.lev.marq <- lev.marq()
solve.nelder.mead <- nelder.mead()
solve.sim.anneal <- sim.anneal()

solve.bfgs.deltas <- reshape.with.diag(solve.bfgs$par, n_v)
solve.lev.marq.deltas <- reshape.with.diag(solve.lev.marq$par, n_v)
solve.nelder.mead.deltas <- reshape.with.diag(solve.nelder.mead$par, n_v)
solve.sim.anneal.deltas <- reshape.with.diag(solve.sim.anneal$par, n_v)

signif(solve.deltas, digits=3)

#############################################
## Save workspace. Can update markdown here.
#############################################
save.image("full_simulation.RData")

update.readme = function()
{
  library(knitr)
  setwd("../doc")
  knit2html("README.Rmd")
  setwd("../bin")
}