# Biplot plotting function

Creates a multiple correspondence analysis (MCA) biplot `biplFig()` was
renamed to [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Usage

``` r
# S3 method for class 'missmi'
plot(
  x,
  engine = c("ggplot2", "base"),
  Z.col = "#61223b",
  CLP.col = "#b79962",
  Z.pch = 19,
  CLP.pch = 15,
  Z.cex = 1.5,
  CLP.cex = 1.7,
  title = NULL,
  ...
)

biplFig(
  x,
  engine = c("ggplot2", "base"),
  Z.col = "#61223b",
  CLP.col = "#b79962",
  Z.pch = 19,
  CLP.pch = 15,
  Z.cex = 1.5,
  CLP.cex = 1.7,
  title = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `missmi`, obtained from preceding function
  [`missmi()`](https://jnienk.github.io/GPAbin/reference/missmi.md).

- engine:

  Select whether the plot should be generated with `ggplot2` or `base`,
  with default `ggplot2`.

- Z.col:

  Colour of sample coordinates

- CLP.col:

  Colour of category level point coordinates

- Z.pch:

  Plotting character of sample coordinates

- CLP.pch:

  Plotting character of category level point coordinates

- Z.cex:

  Size of plotting character for sample points

- CLP.cex:

  Size of plotting character for category level point points

- title:

  Title of the GPAbin biplot.

- ...:

  additional arguments.

## Value

- If `compdat = NULL` in
  [`evaluation`](https://jnienk.github.io/GPAbin/reference/evaluation.md),
  only a GPAbin biplot will be constructed.

- If a complete data set (`compdat`) was specified in
  [`evaluation`](https://jnienk.github.io/GPAbin/reference/evaluation.md),
  two biplots will be constructed: (1) Complete MCA biplot and (2)
  GPAbin biplot.

## Examples

``` r
data(implist)
missbp <- missmi(implist)|> DRT() |> GPAbin() |> plot()


# \donttest{
data(implist)
data(compdat)
missbp <- missmi(implist)|> DRT() |> GPAbin() |> 
evaluation(compdat = compdat) |> plot()



### GPAbin biplot

missbp$plot


### MCA biplot

missbp$plotC# }

```
