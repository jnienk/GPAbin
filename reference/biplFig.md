# Biplot function

Creates a multiple correspondence analysis (MCA) biplot

## Usage

``` r
biplFig(
  missbp,
  Z.col = "#61223b",
  CLP.col = "#b79962",
  Z.pch = 19,
  CLP.pch = 15,
  Z.cex = 1.5,
  CLP.cex = 1.7,
  title = ""
)
```

## Arguments

- missbp:

  An object of class `missbp` obtained from preceding function
  [`missmi()`](https://jnienk.github.io/GPAbin/reference/missmi.md)

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

  Title of the plot

## Value

- If `compdat = NULL` in
  [`evalMeas`](https://jnienk.github.io/GPAbin/reference/evalMeas.md),
  only a GPAbin biplot will be constructed.

- If a complete data set (`compdat`) was specified in
  [`evalMeas`](https://jnienk.github.io/GPAbin/reference/evalMeas.md),
  two biplots will be constructed: (1) Complete MCA biplot and (2)
  GPAbin biplot.

## Examples

``` r
data(implist)
missbp <- missmi(implist)|> DRT() |> GPAbin() |> biplFig()
```
