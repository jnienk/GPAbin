# Generic print function for objects of class missmi

This function is used to print output when the missmi biplot object is
created.

## Usage

``` r
# S3 method for class 'missmi'
print(x, ...)
```

## Arguments

- x:

  an object of class `missmi`.

- ...:

  additional arguments.

## Value

This function will not produce a return value, it is called for side
effects.

## Examples

``` r
data(missdat)
missbp <- missmi(missdat)
data(implist)
missbp <- missmi(implist)
print(missbp)
#> [1] "There are 5 imputations / variations of your data available."
```
