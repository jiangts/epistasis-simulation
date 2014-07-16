### Looking at the little cross dataset
#### Methods:

1. Used `cape` method `read.population` to load the dataset

2. Chose variants to analyze based off of a `cape` analysis. Ran the analysis, but with a small number of permutations (only 3)  

3. Plotted network, and found underlying data for network in `data.obj$var.to.var.p.val`
```
##      Source Target     Effect        SE |Effect|/SE P_empirical p.adjusted
## [1,]     87     97 -1.1823088 0.1945755    6.076350           0          0   
## [2,]     25     22 -0.7936384 0.1406786    5.641502           0          0   
## [3,]     89     45 -0.9320404 0.1675941    5.561296           0          0   
## [4,]     42     45 -0.8505218 0.1616167    5.262587           0          0   
## [5,]      5     25 -0.9196543 0.1811451    5.076892           0          0   
## [6,]     68     90 -0.7551244 0.1548268    4.877221           0          0   
```

4. Chose all interactions for which adjusted p value is 0, and plotted the resulting variant network in python using `networkx`.

5. Found subcomponents of the resulting graph to analyze with the multi-variant cape method.
```
>>> import net 
>>> G=net.makeNetwork()
>>> net.connectivity(G)
how many connected components? nx.number_connected_components(UG)
8
what are the subgraphs? nx.connected_components(UG)
[[0, 1, 2, 4, 11, 12, 21, 24, 33, 36, 42, 43, 45, 50, 53, 55, 58, 61, 62, 66, 68, 69, 70, 86, 87, 88, 89, 90, 97, 98, 99], [25, 22, 5, 6, 7], [75, 18, 83, 76, 23], [3, 37], [35, 14], [19, 92], [80, 28], [54, 30]]
```

6. Decided to analyze the 5 variant network consisting of markers 75, 18, 83, 76 and 23, and as a test, analyze the 2 variant network with markers 3 and 37.

5 variant analysis: 
```
##            [,1]       [,2]        [,3]       [,4]       [,5]
## [1,]  0.0000000 -0.1083341 -0.11851270  0.4454975  0.1229988
## [2,]  1.1850168  0.0000000  0.51538936 -0.7563358 -4.2744015
## [3,]  1.6747963 -1.4770575  0.00000000 -1.4563781  3.1019369
## [4,]  1.6511558  1.3575885  0.04438425  0.0000000  0.3988067
## [5,] -0.3802949  0.5437068  0.07237885  0.3202711  0.0000000
```

This does indeed put edges in all the appropriate locations and directions, and additionally infers a host of additional interactions (but without significance tests). However, it includes positive interaction values.

2 variant analysis:
```
## > solve.bfgs.deltas
##            [,1]       [,2]
## [1,]  0.0000000 -0.1593597
## [2,] -0.6494814  0.0000000
```
```
## > solve.nelder.mead.deltas
##          [,1]     [,2]
## [1,] 0.000000 -3.24637
## [2,] 1.040651  0.00000
```

Problems:
* BFGS and Levenberg Marquardt are no longer converging to the same optimum. Furthermore, the function evaluation is around 2, whereas a plain guess of .5's for all the parameter values gives a value around 2.5. In the simulation case, the function evaluation was nearly 0.
* There are no positive interaction values in the original `cape` analysis...?
* Concerningly, the local minima seems to be pretty powerful even in the 2 variant case...
