# Function to unify coordinates of multiple configurations

Combines multiple configurations from dimension reduction solutions
applied to multiple imputed data sets

## Usage

``` r
GPAbin(missbp, G.target = NULL)
```

## Arguments

- missbp:

  An object of class `missmi` obtained from preceding function
  [`missmi()`](https://jnienk.github.io/GPAbin/reference/missmi.md)

- G.target:

  Target configuration. If not specified the centroid configuration will
  be used as the target.

## Value

The `missbp` object is appended with the following objects:

- Z.GPA.list:

  List containing the sample coordinates for each MI after GPA

- CLP.GPA.list:

  List containing the CLPs for each MI after GPA

- G.target:

  Target configuration

- Z.GPAbin:

  Sample coordinates for the GPAbin biplot

- CLP.GPAbin:

  CLPs for the GPAbin biplot

See also
[`missmi`](https://jnienk.github.io/GPAbin/reference/missmi.md),
[`impute`](https://jnienk.github.io/GPAbin/reference/impute.md) and
[`DRT`](https://jnienk.github.io/GPAbin/reference/DRT.md).

For more detail, refer to Nienkemper-Swanepoel, J., le Roux, N. J., &
Gardner-Lubbe, S. (2021). GPAbin: unifying visualizations of multiple
imputations for missing values. Communications in Statistics -
Simulation and Computation, 52(6), 2666–2685.
https://doi.org/10.1080/03610918.2021.1914089.

## Examples

``` r
data(implist)
missbp <- missmi(implist) |> DRT() |> GPAbin()
```
