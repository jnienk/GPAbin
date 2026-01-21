# Dimension reduction

Multiple correspondence analysis is performed on the multiple imputed
datasets

## Usage

``` r
DRT(missbp, method = c("MCA"))
```

## Arguments

- missbp:

  An object of class `missbp` obtained from preceding function
  [`missmi()`](https://jnienk.github.io/GPAbin/reference/missmi.md)

- method:

  Select a dimension reduction technique. In the current version `MCA`
  is available.

## Value

The `missbp` object is appended with the following objects:

- Z:

  List of sample coordinates

- CLP:

  List of category level point coordinates

- lvls:

  List of category level names

- m:

  Number of multiple imputations

See also [`missmi`](https://jnienk.github.io/GPAbin/reference/missmi.md)
and [`impute`](https://jnienk.github.io/GPAbin/reference/impute.md).

## Examples

``` r
data(implist)
missbp <- missmi(implist) |> DRT()
```
