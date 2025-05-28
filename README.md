# **GPAbin**
GPAbin biplots to unify visualisations from multiple imputed data.

## Getting started

To install the package from Github and load using

```
remotes::install_github("jnienk/GPAbin")
library(GPAbin)
```

## About the name

Generalised orthogonal Procrustes analysis (**GPA**) is used to align individual biplots before combining their separate coordinate sets into an average coordinate matrix, to mimic Ru**bin**'s rules. 
Finally, this average coordinate matrix is then utilized to construct a single biplot called a GPAbin biplot. 