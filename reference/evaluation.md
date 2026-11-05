# Evaluation measures when complete data is available

Calculates measures of comparison based on distances between two
configurations in two dimensions.

## Usage

``` r
evaluation(missbp, compdat = NULL)

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

- tab:

  Returns a table with five evaluation measures: 1-Procrustes Statistic
  (1-PS), Similarity Proportion (SP), Response Profile Recovery (RPR),
  Absolute Mean Bias (AMB), Root Mean Squared Bias (RMSB)

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

**Interpretation** GPAbin preserves the distances of the original biplot
when large values (close to 1) are obtained for 1-PS, SP and RPR and
small values for AMB and RMSB relative to other bias measures.

For more detail, refer to Nienkemper-Swanepoel, J., le Roux, N. J., &
Gardner-Lubbe, S. (2021). GPAbin: unifying visualizations of multiple
imputations for missing values. Communications in Statistics -
Simulation and Computation, 52(6), 2666–2685.
https://doi.org/10.1080/03610918.2021.1914089.

## See also

[`missmi`](https://jnienk.github.io/GPAbin/reference/missmi.md),
[`impute`](https://jnienk.github.io/GPAbin/reference/impute.md),
[`DRT`](https://jnienk.github.io/GPAbin/reference/DRT.md),
[`GPAbin`](https://jnienk.github.io/GPAbin/reference/GPAbin.md)

## Examples

``` r
# \donttest{
data(compdat)
data(implist)
missbp <- missmi(implist) |> DRT() |> GPAbin() |> evalMeas(compdat=compdat)
missbp$tab# }
#> 
#> 
#> |     | Evaluation measures|
#> |:----|-------------------:|
#> |1-PS |              0.9873|
#> |SP   |              0.9646|
#> |RPR  |              0.8720|
#> |AMB  |              0.1268|
#> |RMSB |              0.1554|
```
