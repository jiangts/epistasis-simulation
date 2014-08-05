setwd("~/JAX/simulation/v1.5/bin/lab/")
library("MASS")
source("helper_fn.R")
load("steam_data.RData")

X<-little.cross$geno.for.pairscan[,-102] #* 2
Y<-little.cross$ET

covars = 101 ##Sex Covar
#,c(25, 22, 5, 6, 7)
#,c(18,23,76,83,75)
#,c(3,37)
#,c(32,33,50,51,52)
data = get.relevant.markers(c(6,21,59,87,100, covars)) ##Sex Covar
X = data$X
Y = data$Y

sex = X[,as.character(covars)] ##super janky. Grabbing sex off of X ##Sex Covar
X = X[,-(ncol(X))] # Then deleting it (because of get.design.mat) ##Sex Covar

# need to define n_v
n_v <- ncol(X)

design = get.design.mat(X,n_v)
design = cbind(design, sex) # Then re-adding it ##Sex Covar

solve.betas <- ginv(design) %*% Y
solve.betas = solve.betas[-(nrow(solve.betas)),] ##Sex Covar
solve.all.betas = solve.betas

guess = rep(.2, n_v^2-n_v)

solve.lev.marq <- lev.marq(guess)

solve.m.vals <- infer.m(solve.lev.marq$par)

#with covariate
#            [,1]      [,2]       [,3]
#[1,]  0.00000000 -1.085573 -0.8469582
#[2,]  0.04155386  0.000000  0.2367756
#[3,] -0.53430189 -2.590950  0.0000000 cost: 0.0003856
#no covariate
#           [,1]       [,2]      [,3]
#[1,]  0.0000000 -0.9004915 -1.209559
#[2,] -0.3059209  0.0000000 -1.106612
#[3,]  0.4784661 -0.1887149  0.000000  cost: 0.0007612

distrib.size = 5000

#CV.dist = generate.CV.distribution(design, maxerr=.25, cv.sample.size = distrib.size, cv.sample.frac=.8)
perm.dist = generate.perm.distribution(design, maxerr=70, perm.sample.size = distrib.size)

#p.values = get.p.values(CV.dist)
p.values = get.p.values(perm.dist, solve.m.vals$par)
#TODO: instead of assuming a t distribution, get percentage of values past 0. 
#It's pretty cheap to get examples.
#number.bad = get.bad.occurrences(CV.dist)
number.bad = get.bad.perm.occurrences(perm.dist, solve.m.vals$par)
p.vals = number.bad / distrib.size

print(p.vals)
for(i in 1:ncol(CV.dist))
{
  hist(CV.dist[,i], breaks = 20)
}

data = CV.dist
data <- apply(data, 2, scale)

a = apply(CV.dist, 2, function(x){(0 - mean(x)) / sd(x)})
s <- apply(data, 2, sd)
n <- nrow(data)
xbar <- colMeans(data)
t <- (xbar-a)/(s/sqrt(n))
p.values = 2*pt(-abs(t),df=n-1)

## Manually scaling
#(x - mean(x)) / sd(x)

