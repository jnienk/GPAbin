# Category level prediction

Predicts category levels from an MCA based biplot using the distances
between coordinates

## Usage

``` r
CLPpred(CLPs = CLPs, Zs = Zs, p = p, n = n, lvls = lvls, datIN = datIN)
```

## Arguments

- CLPs:

  Category level point coordinates

- Zs:

  Sample coordinates

- p:

  Number of variables

- n:

  Number of samples

- lvls:

  Names of category levels

- datIN:

  Input data from which `CLPs` and `Zs` are obtained

## Value

- predCL:

  Final predicted categorical data set
