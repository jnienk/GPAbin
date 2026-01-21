# First step before constructing unified biplots

This function produces a list of elements to be used when producing a
GPAbin biplot.

## Usage

``` r
missmi(data)
```

## Arguments

- data:

  input data frame or list

## Value

- X:

  The processed data

- m:

  Number of multiple imputations applied

- n:

  The number of samples

- p:

  The number of variables

- miss_pct:

  Percentage of missing values

## Examples

``` r
data(missdat)
missbp <- missmi(missdat)
data(implist)
missbp <- missmi(implist)
```
