################ Helper Functions ################

##################################################
# prec: list of markers to subset
# postc: filtered marker data with all NAs removed
##################################################
get.relevant.markers = function(l)
{
  X.anal <- X[,l]
  idxs = which(!is.na(rowSums(X.anal)))
  X.anal <- X.anal[idxs,]
  Y.anal <- Y[idxs,]
  return(list("X"=X.anal,"Y"=Y.anal))
}

##################################################
# prec: input raw selected marker data
# postc: return design matrix with interaction term columns
##################################################
get.design.mat = function(X.anal,n_v)
{
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
  return(cbind(1, X.anal))
}

##################################################
# prec: an optimization result
# postc: a pretty form of the delta values matrix
##################################################
pretty.deltas = function(res)
{
  print(reshape.with.diag(res$par, n_v))
}

##############################
##############################
##################################################
# prec: create a sample CV population and calculate beta values
# postc: TEMPORARILY mutate solve.betas. Horrible design, but I 
#        don't know how to do it better.
##################################################
generate.CV.sample.old = function(design, cv.sample.frac)
{
  cv.size = ceiling(nrow(design)*cv.sample.frac)
  cv.sample = sample(1:nrow(design), cv.size, replace = FALSE)
  s.X = design[cv.sample,]
  x.Y = Y[cv.sample,]
  solve.betas <<- ginv(s.X) %*% x.Y
  #print(tail(cv.sample))
} ##############################
##############################
##############################

##################################################
# prec: create a sample CV population and calculate beta values
# postc: TEMPORARILY mutate solve.betas. Horrible design, but I 
#        don't know how to do it better.
##################################################
generate.CV.sample = function(design, cv.sample.frac)
{
  cv.size = ceiling(nrow(design)*cv.sample.frac)
  cv.sample = sample(1:nrow(design), cv.size, replace = FALSE)
  s.X = design[cv.sample,]
  x.Y = Y[cv.sample,]
  beta.sex.covar = ginv(s.X) %*% x.Y
  solve.betas <<- beta.sex.covar[-(nrow(beta.sex.covar)),]
}

##################################################
# prec: Input design matrix, sample size, CV sample fraction size, and maxerr
# postc: return distribution of resulting delta values.
##################################################
generate.CV.distribution = function(design, cv.sample.size = 50, cv.sample.frac=.6, 
                                    maxerr = .001)
{
  CV.dist = matrix(NA, 1, n_v^2-n_v)
  for(i in 1:cv.sample.size)
  {
    generate.CV.sample(design, cv.sample.frac)
    #res <- lev.marq()
    res <- infer.m(lev.marq()$par)
    CV.dist = rbind(CV.dist, res$par)
    if(res$deviance > maxerr)
    {
      cat("Warning: optimization err greater than maxerr. Sample Nu. ", i)
    }
  }
  solve.betas <<- solve.all.betas
  return(CV.dist[-1,])
}


#######################################################
#######################################################
# Screw it, I'm doing permutation tests!

##################################################
# prec: data
# postc: permute each column of data separately
##################################################
permuteData = function(data) {
  for(i in 1:ncol(data)){
    data[,i] = data[sample(nrow(data)),i]
  }
  return(data)
}

##################################################
# prec: create a sample CV population and calculate beta values
# postc: TEMPORARILY mutate solve.betas. Horrible design, but I 
#        don't know how to do it better.
##################################################
generate.perm.sample = function(design)
{
  #don't permute the sex covar with it! ugh
  sex.orig.ordering = design[,"sex"]
  s.X = permuteData(design)
  #s.X = design[sample(nrow(design)),]
  s.X[,"sex"]=sex.orig.ordering
  beta.sex.covar = ginv(s.X) %*% Y
  solve.betas <<- beta.sex.covar[-(nrow(beta.sex.covar)),]
}

##################################################
# prec: Input design matrix, sample size, CV sample fraction size, and maxerr
# postc: return distribution of resulting delta values.
##################################################
generate.perm.distribution = function(design, perm.sample.size = 50, 
                                    maxerr = .001)
{
  perm.dist = matrix(NA, 1, n_v^2-n_v)
  for(i in 1:perm.sample.size)
  {
    generate.perm.sample(design)
    #res <- lev.marq()
    res <- infer.m(lev.marq()$par)
    if(res$deviance > maxerr)
    {
      cat("\nWarning: ", res$deviance, 
          " = opt err > maxerr = ", maxerr, ". Sample Nu. ", i)
    } else {
      perm.dist = rbind(perm.dist, res$par)
    }
  }
  solve.betas <<- solve.all.betas
  return(perm.dist[-1,])
}


#######################################################
#######################################################







##################################################
# prec: input all CV distribution of delta values
# postc: return p values of each distribution against 0
##################################################
get.p.values = function(data, a = rep(0, ncol(data)))
{
  #a <- solve.lev.marq$par
  s <- apply(data, 2, sd)
  n <- nrow(data)
  xbar <- colMeans(data)
  t <- (xbar-a)/(s/sqrt(n))
  cat("t values are: ", t)
  return(2*pt(-abs(t),df=n-1))
}

##################################################
# prec: input all CV distribution of delta values
# postc: return "bad occurrences" of data
##################################################
# a "bad occurrence" is a sample that goes against the story told by our distribution.
# so, if our distribution says them mean is 2, but 1/3 of the samples are <=0, we have lots of
# bad occurrences.
get.bad.occurrences = function(data, ref = 0)
{
  xbar = colMeans(data)
  gtRef = colSums(data >= ref)
  ltRef = colSums(data <= ref)
  out = matrix(NA, 1, ncol(data))
  for(i in 1:length(xbar))
  {
    if(xbar[i] > ref){
      out[i] = ltRef[i]
    }
    if(xbar[i] < ref){
      out[i] = gtRef[i]
    }
  }
  return(out)
}

##################################################
# prec: input all CV distribution of delta values
# postc: return "bad occurrences" of data
##################################################
# a "bad occurrence" is a sample that goes against the story told by our distribution.
# so, if our distribution says them mean is 2, but 1/3 of the samples are <=0, we have lots of
# bad occurrences.
get.bad.perm.occurrences = function(data, mean)
{
  out = matrix(NA, 1, ncol(data))
  xbar = colMeans(data)
  for(i in 1:length(mean))
  {
    gtRef = colSums(data >= mean[i])
    ltRef = colSums(data <= mean[i])
  
    if(xbar[i] > mean[i]){
      out[i] = ltRef[i]
    }
    if(xbar[i] < mean[i]){
      out[i] = gtRef[i]
    }
  }
  return(out)
}






##################################################
# prec: take arbitrary mat. M and route length k
# postc: return adj. mat. of Hamiltonian paths of length k
##################################################

hamiltonian.len = function(M, k)
{
  if(k-1 == 0) { return(M) }
  out <- M
  for(i in 1:(k-1))
  {
    out <- out %*% M
    diag(out) <- 0
  }
  return(out)
}

##################################################
# desc: store all combos into a global var 
# purp: this will hopefully save time on optim
##################################################

create.combos <- function(n_v)
{
  out <- list()
  for(k in 1:n_v)
  {
    out[[k]] <- t(combn(c(1:n_v), k))
  }
  return(out)
}

##################################################
# prec: take interaction matrix A 
# postc: return interaction matrix I
##################################################

get.I.from.A = function(A)
{
  n_v <- nrow(A)
  I <- diag(n_v+1)
  
  if(!exists("all.combos") || length(all.combos) != n_v){
    all.combos <<- create.combos(n_v)
  }
  for(k in 2:n_v)
  {
    combos <- all.combos[[k]]
    deltas <- matrix(0, nrow(combos), n_v)
    for(i in 1:nrow(combos))
    {
      combo <- combos[i,]
      active_vars <- A[combo,combo]
      routes <- hamiltonian.len(active_vars, k-1)
      delta_row <- colSums(routes) 
      deltas[i,combo] <- delta_row
    }
    I <- rbind(I, cbind(matrix(0,nrow(deltas),1), deltas))
  }
  return(I)
}

##################################################
# prec: take interaction matrix A and main effects
# postc: return full beta matrix
##################################################

get.betas.from.network = function(A, ME) 
{
  I <- get.I.from.A(A)
  return(I %*% ME)
}

### get_interactions.R's fns

##################################################
# prec: take unrolled delta guess vector
# postc: return squared error cost
##################################################

#TODO: poor coding; grabbing n_v and solve.betas from global scope
cost_fn <- function(guess)
{ 
  #shape into matrix and insert diagonal of 0's
  A <- reshape.with.diag(guess, n_v)
  betas.guess <- get.betas.from.network(A, solve.betas[1:(n_v+1),]) #nice code reusage!...
  
  return(sum((solve.betas - betas.guess)^2))
}

library("minpack.lm")
lma_cost_fn <- function(guess)
{ 
  #shape into matrix and insert diagonal of 0's
  A <- reshape.with.diag(guess, n_v)
  betas.guess <- get.betas.from.network(A, solve.betas[1:(n_v+1),]) #nice code reusage!...
  
  return(solve.betas - betas.guess)
}

##################################################
# prec: take unrolled vector and side length
# postc: return re-rolled matrix with a diagonal of 0's
##################################################

reshape.with.diag <- function(unr, s)
{
  ind <- s*0:(s-1)
  val <- c(unr, rep(0, length(ind)))
  id <- c(seq_along(unr), ind+0.5)
  return(matrix(val[order(id)], s, s))
}

#############################################
# TODO: calculate gradient as well...
# Optimization methods
#############################################

bfgs = function(par = rep(.5, n_v^2-n_v))
{
  res <- optim(par, cost_fn, NULL, method = "BFGS")
  return(res)
}

nelder.mead = function(par = (1:(n_v^2-n_v) * .1))
{
  res <- optim(par, cost_fn, method = "Nelder-Mead",
               control = list(maxit = 20000))
  return(res)
}

sim.anneal = function(par = rep(.5, n_v^2-n_v))
{
  res <- optim(par, cost_fn, method = "SANN",
               control = list(maxit = 20000, temp = 20, parscale = rep(.5, n_v^2-n_v) ))
  return(res)
}

lev.marq = function(par = rep(.5, n_v^2-n_v))
{
  nls.out <- nls.lm(par, lower=NULL, upper=NULL, lma_cost_fn, jac = NULL,
                control = nls.lm.control())
  return(nls.out)
}

################ Helper Functions ################

##################################################
# prec: take optimized delta values
# postc: infer interaction values m
##################################################
infer.m = function(delta_result, par = rep(.5, n_v^2-n_v)){
  I = get.I.from.A(reshape.with.diag(delta_result, n_v)) 
  Deltas <<- I[0:-n_v-1,-1] #isolate the delta sums
  
  flat.combos <<- flatten.combos(all.combos)
  nls.out <- nls.lm(par, lower=NULL, upper=NULL, delta.to.m.cost.fn.lma, jac = NULL,
                    control = nls.lm.control())
  #m.res <- optim(rep(.5, n_v^2-n_v), delta.to.m.cost.fn, NULL, method = "BFGS")
  return(nls.out)
}

flatten.combos = function(ac)
{
  ac = ac[-1] #remove 1 variant active cases.
  out = list()
  for(i in 1:length(ac))
  {
    for(j in 1:nrow(ac[[i]]))
    {
      out = c(out, list(ac[[i]][j,]))
    }
  }
  return(out)
}

is.subset <- function(x, y) all(x %in% y)
reshape.with.neg.one.diag <- function(unr, s) #super tricky!
{
  ind <- s*0:(s-1)
  val <- c(unr, rep(-1, length(ind)))
  id <- c(seq_along(unr), ind+0.5)
  return(matrix(val[order(id)], s, s))
}
delta.to.m.cost.fn = function(M){  #OMIT ARGS Deltas, flat.combos. Again, crap coding style
  M = reshape.with.neg.one.diag(M, n_v)
  
  total.sq.err = 0
  for(i in 1:length(flat.combos))
  {
    ids <- numeric()
    for(j in 1:i) #see if previous combos are subsets and include them.
    {
      if(is.subset(flat.combos[[j]], flat.combos[[i]]) == TRUE)
      {
        ids = c(ids, j)
      }
    }
    combo = flat.combos[[i]]
    active.deltas = colSums( rbind(1, Deltas[ids, combo]) )
    active.m = M[combo, combo]
    #the tricky step!
    all.errors = t(matrix(active.deltas)) %*% active.m + 1
    total.sq.err = total.sq.err + sum(all.errors^2)
  }
  return(total.sq.err)
}


delta.to.m.cost.fn.lma = function(M){  #OMIT ARGS Deltas, flat.combos. Again, crap coding style
  M = reshape.with.neg.one.diag(M, n_v)
  
  all.err.values = c()
  for(i in 1:length(flat.combos))
  {
    ids <- numeric()
    for(j in 1:i) #see if previous combos are subsets and include them.
    {
      if(is.subset(flat.combos[[j]], flat.combos[[i]]) == TRUE)
      {
        ids = c(ids, j)
      }
    }
    combo = flat.combos[[i]]
    active.deltas = colSums( rbind(1, Deltas[ids, combo]) )
    active.m = M[combo, combo]
    #the tricky step!
    all.errors = t(matrix(active.deltas)) %*% active.m + 1
    all.err.values = c(all.err.values, all.errors)
  }
  #print(length(all.err.values))
  return(all.err.values)
}


