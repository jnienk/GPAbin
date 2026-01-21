# Generalised Orthogonal Procrustes Analysis

This function contains the OPA function to compare two configurations
and the GPA function for multiple configuration comparisons

## Usage

``` r
GPA(Xk, G.target = NULL, iter = 500, eps = 0.001)
```

## Arguments

- Xk:

  list containing the testee configurations which is updated on \#each
  iteration

- G.target:

  Target configuration. If not specified the centroid configuration will
  be used as the target

- iter:

  Number of iterations allowed before convergence

- eps:

  Threshold value for convergence of the alogrithm

## Value

- Xk.F:

  List containing the updated testee configurations

- sk.F:

  Vector containing the final scaling factors

- Qk.F:

  List containing the final rotation matrices

- Gmat:

  Final target configuration

- sum.sq:

  Final minimised sum of squared distance
