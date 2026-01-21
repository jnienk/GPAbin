# Multiple imputation

Choose between four available multiple imputation strategies in `R`.

## Usage

``` r
impute(missbp, imp.method = c("MIMCA", "jomo", "DPMPM", "mice"), m = 5)
```

## Arguments

- missbp:

  An object of class `missmi` obtained from preceding function
  [`missmi()`](https://jnienk.github.io/GPAbin/reference/missmi.md).

- imp.method:

  Select one of four imputation methods: `MIMCA`, `jomo`, `DPMPM`,
  `mice`

- m:

  Number of multiple imputations

## Value

The `missbp` object is appended with the following object:

- dataimp:

  List of imputed data

See also [`MIMCA`](https://rdrr.io/pkg/missMDA/man/MIMCA.html),
[`jomo1cat`](https://rdrr.io/pkg/jomo/man/jomo1cat.html),
[`mi`](https://rdrr.io/pkg/mi/man/04mi.html) and
[`mice`](https://amices.org/mice/reference/mice.html).

## Examples

``` r
# \donttest{
data(missdat)
missbp <- missmi(missdat) |> impute(imp.method="DPMPM", m=5)# }
```
