################ Helper Functions ################

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

##################################################
# prec: take optimized delta values
# postc: infer interaction values m
##################################################

infer.m = function(D){
  I = get.I.from.A(reshape.with.diag(D, n_v))
  
  #total variant activity
  As = colSums(I)[-1] #column sums and remove first entry
}