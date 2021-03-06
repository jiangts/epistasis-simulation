Extending R/cape: Multivariant Analysis
========================================================

The purpose of this project is to test a new model to extend `R/cape` for the analysis of more than 2 variants simulataneously.

There are two parts to the code:
* `create_data.R`, where we fabricate a dataset given an underlying variant-to-variant interaction network, and
* `get_interactions.R` where we use our new model on the fabricated marker and phenotype data to re-infer the underlying variant interactions.

Using this methodology, we can implement the new model of analysis while testing its robustness as it responds to varying levels of noise that we can systematically introduce into the data.

------

Running the code: an example
========================================================
## create_data.R
Here is an example of the data generated by the script `create_data.R`. We set several parameters in the fabrication of this dataset, most notably
* 1000 individuals,
* 3 phenotypes,
* 3 variants,
* 50% allelic frequency, and
* __No noise__ -- this will be updated shortly

Then, using the underlying interaction matrix $\boldsymbol{A}$, set to
<!-- 

```r
library("MASS")
library("Matrix")
load("~/JAX/simulation/v1.2/bin/full_simulation.RData")
```

```
## Warning: cannot open compressed file
## '/Users/mingya/JAX/simulation/v1.2/bin/full_simulation.RData', probable
## reason 'No such file or directory'
```

```
## Error: cannot open the connection
```
 -->

```r
print(A)
```

```
##      [,1] [,2] [,3]
## [1,]  0.0  0.4  0.0
## [2,]  0.1  0.0  0.3
## [3,]  0.0  0.0  0.0
```

and the main effect matrix

```r
print(ME_betas)
```

```
##      [,1] [,2] [,3]
## [1,]  3.0  3.5  4.0
## [2,]  3.5  4.0  3.0
## [3,]  4.0  3.0  3.5
```

we created the $\beta$ matrix.

```r
print(betas)
```

```
##      [,1] [,2] [,3]
## [1,] 1.00 1.00 1.00
## [2,] 3.00 3.50 4.00
## [3,] 3.50 4.00 3.00
## [4,] 4.00 3.00 3.50
## [5,] 1.70 1.95 1.60
## [6,] 0.00 0.00 0.00
## [7,] 1.20 0.90 1.05
## [8,] 0.48 0.36 0.42
```

Now, using the linear regression model
$$
Y \sim X \cdot \beta + \epsilon,
$$
where $Y$ is the phenotype matrix (columns represent phenotypes, rows invididuals), $X$ is the genotype matrix, and $\beta$ is the weights matrix, we get the following fabricated data:


```r
colnames(X) <- c("x0", "x1", "x2", "x3", "x12", "x13", "x23", "x123")
head(X)
```

```
##      x0 x1 x2 x3 x12 x13 x23 x123
## [1,]  1  2  2  2   4   4   4    8
## [2,]  1  2  1  0   2   0   0    0
## [3,]  1  1  1  1   1   1   1    1
## [4,]  1  0  1  2   0   0   2    0
## [5,]  1  2  2  0   4   0   0    0
## [6,]  1  0  1  1   0   0   1    0
```

```r
head(Y)
```

```
##       [,1]  [,2]  [,3]
## [1,] 37.44 36.28 35.96
## [2,] 13.90 15.90 15.20
## [3,] 14.88 14.71 14.57
## [4,] 14.90 12.80 13.10
## [5,] 20.80 23.80 21.40
## [6,]  9.70  8.90  8.55
```

------
## get_interactions.R

<!-- 
__TODO: there is currently a problem with delta matrix $\boldsymbol{I}$ generation involving routing and feedback.__ The part in `helper_fn.R` in the `get.I.from.A` method:
```
...
  active_vars <- A[combo,combo]
  routes <- active_vars %^% (k-1) #routes of length k-1
  diag(routes) <- 0
  delta_row <- colSums(routes) 
...
``` -->


To infer the values of the interactions from the marker data $X$ and phenotype data $Y$, we first find the weights $\beta$ of our linear regression model. Since

$$\beta \approx X^{-1} \cdot Y,$$

in `R` we do:

```r
ginv(X) %*% Y
```

```
##           [,1]      [,2]      [,3]
## [1,] 1.000e+00 1.000e+00 1.000e+00
## [2,] 3.000e+00 3.500e+00 4.000e+00
## [3,] 3.500e+00 4.000e+00 3.000e+00
## [4,] 4.000e+00 3.000e+00 3.500e+00
## [5,] 1.700e+00 1.950e+00 1.600e+00
## [6,] 3.905e-14 3.987e-14 4.059e-14
## [7,] 1.200e+00 9.000e-01 1.050e+00
## [8,] 4.800e-01 3.600e-01 4.200e-01
```

Notice that these values are essentially equal to the original $\beta$ values. That's because we didn't add any noise to the data yet.

Continuing, we want to find the interaction matrix $\boldsymbol{A}$ filled with $\delta$'s that best fits the following system by the least squares criterion.

$$\begin{cases} 
  \beta_{12}^1 = \beta_{1}^1\delta_{21} + \beta_{2}^1\delta_{12} \\\ 
  \beta_{12}^2 = \beta_{1}^2\delta_{21} + \beta_{2}^2\delta_{12} \\\ 
  \ldots \\\ 
  \beta_{123}^2 = \beta_{1}^2(\delta_{21} + \delta_{31} + \delta_{23}\delta_{31} + \delta_{32}\delta_{21}) + \beta_{2}^2(\delta_{12} + \delta_{32} + \delta_{13}\delta_{32} + \delta_{31}\delta_{12}) \\\ 
  \beta_{123}^3 = \beta_{1}^3(\delta_{21} + \delta_{31} + \delta_{23}\delta_{31} + \delta_{32}\delta_{21}) + \beta_{2}^3(\delta_{12} + \delta_{32} + \delta_{13}\delta_{32} + \delta_{31}\delta_{12}) 
\end{cases}$$

Enumerating all the cases, we see that we have 12 equations and 6 variables. (In general, we will have $p (2^v - v - 1)$ equations and ${v \choose 2}$ variables.)

Now, our task is to find the best $\delta$'s to fit the 12 equations. We use a series of three optimization algorithms to look for the minima: BFGS, a gradient search method, Nelder-Mead, a simplex method, and Simulated Annealing, an entropy-based random global search method.

We show the results for each here:
```
solve.bfgs <- bfgs()
```


```r
print(solve.bfgs)
```

```
## $par
## [1]  1.000e-01  2.632e-10  4.000e-01  3.657e-10 -2.064e-10  3.000e-01
## 
## $value
## [1] 8.085e-19
## 
## $counts
## function gradient 
##       86       40 
## 
## $convergence
## [1] 0
## 
## $message
## NULL
```

```r
signif(reshape.with.diag(solve.bfgs$par, n_v), digits=4)
```

```
##           [,1]      [,2]       [,3]
## [1,] 0.000e+00 4.000e-01 -2.064e-10
## [2,] 1.000e-01 0.000e+00  3.000e-01
## [3,] 2.632e-10 3.657e-10  0.000e+00
```

```
solve.nelder.mead <- nelder.mead()
```


```r
print(solve.nelder.mead)
```

```
## $par
## [1]  9.982e-02  8.360e-05  4.001e-01 -6.436e-05 -3.736e-05  3.000e-01
## 
## $value
## [1] 6.142e-07
## 
## $counts
## function gradient 
##     1457       NA 
## 
## $convergence
## [1] 0
## 
## $message
## NULL
```

```r
signif(reshape.with.diag(solve.nelder.mead$par, n_v), digits=4)
```

```
##           [,1]       [,2]       [,3]
## [1,] 0.0000000  4.001e-01 -3.736e-05
## [2,] 0.0998200  0.000e+00  3.000e-01
## [3,] 0.0000836 -6.436e-05  0.000e+00
```

```
solve.sim.anneal <- sim.anneal()
```


```r
print(solve.sim.anneal)
```

```
## $par
## [1]  0.1881  0.2623  0.3234  0.1271 -0.2809  0.2055
## 
## $value
## [1] 0.2264
## 
## $counts
## function gradient 
##    20000       NA 
## 
## $convergence
## [1] 0
## 
## $message
## NULL
```

```r
signif(reshape.with.diag(solve.sim.anneal$par, n_v), digits=4)
```

```
##        [,1]   [,2]    [,3]
## [1,] 0.0000 0.3234 -0.2809
## [2,] 0.1881 0.0000  0.2055
## [3,] 0.2623 0.1271  0.0000
```

We see that the results in the solution vector (called `par`, for the parameter vector) is quite similar for BFGS and Nelder Mead. This is reassuring as BFGS is a gradient method (and its performance generally depends on the initial guess), whereas Nelder-Mead is a global search method whose ability to find a minimum is less dependent on the intialization point. Thus, our method can show us with some confidence where the global minimum is.

Furthermore, we see that our method has done well in finding the global minimum in this case, as the output matrices are very close to the original seeded interaction matrix.
