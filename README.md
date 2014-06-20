epistasis-simulation
====================
Creating simulated data to test statistical techniques to detect gene-gene interactions.

### Why this?
The idea is to fabricate a dataset where you know the underlying gene-gene interactions so you can test which of the statistical techniques is best able to pick up the interactions you want it to.

### How do I run it?
The main file is `make_data.R`. The code is commented. Run the whole file.

### Analysis:
Our analysis includes 3 variants, 3 phenotypes, and 1000 subjects.
Our interaction network has the following structure:
```
     [,1] [,2] [,3]
[1,]  0.0  0.9  0.0
[2,]  0.4  0.0  0.8
[3,]  0.0  0.0  0.0
```
This means that delta_12 = .9, delta_21 = .4, and delta_23 = .8

We then do a transformation
```R
int_adj <- int_adj + int_adj %^% 2 + int_adj %^% 3 - (diag(n_v) * (int_adj %^% 2))
```
to turn this underlying interaction network into one that includes indirect influences. We do this because we need a network of overall variant to variant influences in order to create our dataset.
The result is:
```
     [,1] [,2] [,3]
[1,] 0.36 0.90 0.72
[2,] 0.40 0.36 0.80
[3,] 0.00 0.00 0.00
```

We also create our main effects, which are as follows:
```
     [,1] [,2] [,3]
[1,]  1.0  1.0  1.0
[2,]  3.0  3.5  4.0
[3,]  3.5  4.0  3.0
[4,]  4.0  3.0  3.5
```

From this, we create our beta matrix, and separately create a genotype distribution in order to create our whole fabricated dataset.
The objects included in the fabricated dataset should include `design` and `phenos`.

Now is time to "reverse solve" for the interaction values using our statistical techniques.
####cape
The first test is to see how well the current methodology fares.
After performing our "reverse solve" algorithm by doing pairwise regressions on TWO phenotypes at a time, and comparing them to our interaction adjacency matrix from which we created the data, we see that our results (as expected) are not bad.
```
> print(signif(delta_cape2, digits=3))
     v1 v2      d21  d12
[1,]  1  2 4.00e-01 0.90
[2,]  1  3 4.44e-16 0.72
[3,]  2  3 9.99e-15 0.80
> int_adj
     [,1] [,2] [,3]
[1,] 0.36 0.90 0.72
[2,] 0.40 0.36 0.80
[3,] 0.00 0.00 0.00
```
However, when running a "reverse solve" by regressing deltas against ALL phenotypes, we get the following result:
```
> print(signif(B_sep, digits=3))
       [,1]  [,2]   [,3]
[1,]  0.000 3.560  1.340
[2,] -2.670 0.000 -0.429
[3,] -0.815 0.904  0.000
```
Oddly, we get a strong negative signal at [2,1], and weaker negative signals at [3,1], and [2,3].

The interpretation of this is still unclear to me.

### separation of overall influences
In this method, we take all the overall beta outputs given to us by cape, and try to infer the direct influences from those overal influences.
```
> print(signif(A_sep, digits=3))
       [,1]  [,2]    [,3]
[1,]  0.861 0.235  0.2870
[2,] -0.302 0.790  0.3160
[3,]  0.160 0.381 -0.0515
```
Again, this makes little sense to me. For some reason, the [1,1] effect is stronger than the [1,2] effect. The same negative effects remain. This is bad...

### directly search for direct influences using optimization method
In this method, we do the regression on the equation beta_12 = beta1 (d21 + d23 * d31) + beta2 (d12 + d13 * d32) in the 3 variable case. The code does generalize to more variables.
Our output in this case is
```
> dir <- reshape.with.diag(result$par,n_v)
> print(signif(dir, digits=3))
       [,1]  [,2]  [,3]
[1,]  0.000 3.470 0.272
[2,] -2.720 0.000 0.310
[3,]  0.158 0.358 0.000
```
This one still has that strong negative effect on [2,1]. Looks like it never went away.

# Preliminary Conclusion
So far, our ability to detect the interactions using any of the 3 methods are looking pretty bad.
Of course, there are many different ways to turn our underlying interaction matrix into a fabricated dataset, and I hypothesize that the manner by which we create `int_adj` will make a large impact on how well our methods perform. In addition, the statistical techniques themselves can be tweaked.

