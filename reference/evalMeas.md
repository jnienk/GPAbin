# Evaluation measures when complete data is available

Calculates measures of comparison based on distances between two
configurations in two dimensions.

## Usage

``` r
evalMeas(missbp, compdat = NULL)
```

## Arguments

- missbp:

  An object of class `missbp` obtained from preceding function
  [`missmi()`](https://jnienk.github.io/GPAbin/reference/missmi.md).

- compdat:

  Complete data matrix representing the input data of
  [`missmi()`](https://jnienk.github.io/GPAbin/reference/missmi.md)

## Value

The `missbp` object is appended with the following objects:

- eval:

  Returns a data table with five evaluation measures: Procrustes
  Statistic (PS), Similarity Proportion (SP), Response Profile Recovery
  (RPR), Absolute Mean Bias (AMB), Root Mean Squared Bias (RMSB)

- GPApred:

  A dataframe representing predicted categorical responses from the
  GPAbin biplot.

- compPred:

  A dataframe representing predicted categorical responses from the
  complete MCA biplot.

- compZs:

  Sample coordinates for the MCA biplot of the complete data.

- compCLPs:

  CLPs for the MCA biplot of the complete data.

- complvls:

  Names of the CLPs for the MCA biplot of the complete data.

See also
[`missmi`](https://jnienk.github.io/GPAbin/reference/missmi.md),
[`impute`](https://jnienk.github.io/GPAbin/reference/impute.md),
[`DRT`](https://jnienk.github.io/GPAbin/reference/DRT.md) and
[`GPAbin`](https://jnienk.github.io/GPAbin/reference/GPAbin.md).

For more detail, refer to Nienkemper-Swanepoel, J., le Roux, N. J., &
Gardner-Lubbe, S. (2021). GPAbin: unifying visualizations of multiple
imputations for missing values. Communications in Statistics -
Simulation and Computation, 52(6), 2666–2685.
https://doi.org/10.1080/03610918.2021.1914089.

## Examples

``` r
# \donttest{
data(compdat)
data(implist)
missbp <- missmi(implist) |> DRT() |> GPAbin() |> evalMeas(compdat=compdat)# }
```
