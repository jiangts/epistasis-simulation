################ Helper Functions ################

################ 
# desc: turn melted interactions list matrix into adjacency matrix
get.var.adj.matrix = function(interactions, n_v)
{
  #create empty variant interaction adjacency matrix
  int_adj <- matrix(0, nrow = n_v, ncol = n_v)
  
  #populate adjacency matrix
  for (i in 1:nrow(interactions)) {
    int_adj[interactions[i, 1], interactions[i, 2]] = interactions[i, 3]
  }
  return(int_adj)
}
# usage: int_adj <- get.var.adj.matrix(interactions, n_v)

################ 
# desc: create design matrix from number of variants, prob. of "mutant", and n_i
create.genotype.population = function(n_v, p, n_i)
{
  design <- matrix(1, nrow = n_i, ncol = 1)
  
  #get matrix of 0, 1, and 2 using p
  a <- matrix(runif(n_v * n_i, 0, 1), nrow = n_i, ncol = n_v)
  a[a < p] = 1
  a[a != 1] = 0
  b <- matrix(runif(n_v * n_i, 0, 1), nrow = n_i, ncol = n_v)
  b[b < p] = 1
  b[b != 1] = 0
  variants <- a + b
  
  combos <- t(combn(c(1:n_v), 2)) #puts results of n choose 2 in a matrix
  two_variants <- matrix(0, nrow = n_i, ncol = dim(combos)[1])
  #populate two_variants with products.
  for (i in 1:dim(combos)[1])
  { #i is combo iterator
    for(j in 1:n_i)
    { #j is sample iterator
      two_variants[j,i] = variants[j, combos[i,1]] * variants[j, combos[i,2]]
    }
  }
  
  return(cbind(design, variants, two_variants))
}
# usage: design <- create.genotype.population(n_v, var_prob, n_i)

################ 
# desc: get all beta values (need to calculate all interaction beta values)
#want matrix like:
#beta_0^A   beta_0^B
#beta_1^A   beta_1^B
#beta_2^A   beta_2^B
#beta_12^A  beta_12^B
#beta_13^A  beta_13^B
#beta_23^A  beta_23^B
get.beta.values = function(me, vi) 
{ #takes main effects and variant interactions in adj. matrices
  combos <- t(combn(c(1:dim(vi)[1]), 2)) 
  # ^^ do we want to risk running this twice? (in prev fn too)
  
  #currently doing simple model from cape: B_12 = B1*d21 + B2*d12
  beta_ij <- matrix(0, nrow = dim(combos)[1], ncol = dim(me)[2])
  for(i in 1:dim(combos)[1])
  {
    pair <- combos[i,]
    p1 <- pair[1]
    p2 <- pair[2]

    beta_ij[i,pair[1]] <- vi[p2,p1] * me[p1+1,p1] + vi[p1,p2] * me[p2+1,p1]
    beta_ij[i,pair[2]] <- vi[p2,p1] * me[p1+1,p2] + vi[p1,p2] * me[p2+1,p2]
  }
  return(rbind(me, beta_ij))
}
# usage: betas <- get.beta.values(main_effect_betas, int_adj)

################ 
# desc: take betas matrix and number of variants to return matrix B form.
get.beta.adj.matrix = function(betas, n_v) {return(betas[(2+n_v):dim(betas)[1],])}
# usage: B_cape <- get.beta.adj.matrix( ginv(design) %*% phenos )

################ 
# desc: take interaction betas matrix and get cape deltas. output into melted matrix

## THIS IS NOT REGRESSING AGAINST ALL PHENOTYPES. JUST TWO
get.cape.deltas2 = function(betas, n_v)
{
  me <- betas[2:(1+n_v),]
  betas <- get.beta.adj.matrix(betas, n_v)
  combos <- t(combn(1:n_v, 2))
  deltas <- matrix(0,nrow=dim(betas)[1],ncol=2)
  for(i in 1:dim(betas)[1])
  {
    pair <- combos[i,]
    pair_effect <- matrix(c(betas[i,pair[1]], betas[i,pair[2]]))
    weights <- t(me[pair,pair]) #be careful! needs to be transposed!
    deltas[i,] <- t(ginv(weights) %*% pair_effect)
  }
  out <- cbind(combos,deltas)
  colnames(out) <- c("v1", "v2", "d21", "d12")
  return(out)
}
# usage: delta_cape <- get.cape.deltas(B_cape, n_v)

# SAME METHOD THAT REGRESSES DELTA OVER ALL PHENOTYPES
get.cape.deltas = function(betas, n_v)
{
  me <- betas[2:(1+n_v),]
  betas <- get.beta.adj.matrix(betas, n_v)
  combos <- t(combn(1:n_v, 2))
  deltas <- matrix(0,nrow=dim(betas)[1],ncol=2)
  for(i in 1:dim(betas)[1])
  {
    pair <- combos[i,]
    pair_effects <- matrix(betas[i,])
    weights <- t(me[pair,]) #be careful! needs to be transposed!
    deltas[i,] <- t(ginv(weights) %*% pair_effects)
  }
  out <- cbind(combos,deltas)
  colnames(out) <- c("v1", "v2", "d21", "d12")
  return(out)
}
################ 
# desc: convert cape_deltas to adj matrix
get.sep.adj.matrix = function(cape_deltas)
{
  #trying to be cute and generalized... vv
  n_v <- length(union(unique(cape_deltas[,1]),unique(cape_deltas[,2])))
  #create empty variant interaction adjacency matrix
  adj <- matrix(0, nrow = n_v, ncol = n_v)
  
  #populate adjacency matrix
  for (i in 1:nrow(cape_deltas)) {
    adj[cape_deltas[i, 1], cape_deltas[i, 2]] = cape_deltas[i, 4]
    adj[cape_deltas[i, 2], cape_deltas[i, 1]] = cape_deltas[i, 3]
  }
  return(adj)
}
# usage: tmp <- get.sep.adj.matrix(delta_cape)

##################################################
