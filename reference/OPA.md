# Orthogonal Procrustes Analysis

This function performs Orthogonal Procrustes Analysis on centred data

## Usage

``` r
OPA(missbp, compdat, dim = "2D")
```

## Arguments

- missbp:

  An object of class `missmi` obtained from preceding function
  [`missmi()`](https://jnienk.github.io/GPAbin/reference/missmi.md)

- compdat:

  Complete data set, only available for simulated data examples.

- dim:

  Number of dimensions to use in final solutions (`2D` or `All`
  available dimensions.)

## Value

- PS:

  Procrustes Statistic

- RMSB:

  Root Mean Squared Bias

- AMB:

  Absolute Mean Bias

- compZ:

  Sample coordinates representing the complete data set

- compCLP:

  Category level point coordinates representing the complete data set

- complvls:

  Category levels

- compdat:

  Complete data set, only available for simulated data examples
