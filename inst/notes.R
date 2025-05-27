#steps

install.packages(c("devtools", "roxygen2", "testthat", "knitr"))

create_tidy_package(".", copyright_holder = "JNS")
#create_package()

use_mit_license()
#update author name in LICENSE.md

#Rbuild file will ignore the vignette
use_article("GPAbin", "Getting started with GPAbin")

document()

install()
check()
#build_vignettes() #causing problems


library(devtools)
load_all()
document()

#possible work-flow
missmi() |> impute() |> DRT() |> GPAbin() |> plot()
#test

data(missdat)
data(implist)
a <- missmi(missdat)|> impute(imp.method="DPMPM", m=5) |> DRT() |> GPAbin() |> biplFig() |> evalMeas(compdat=comp.sim, dim="2D")
a <- missmi(implist)|> DRT() |> GPAbin() |> biplFig() |> evalMeas(compdat=comp.sim, dim="2D")

#data check
#imputation if needed
#dimension reduction
#unification
#plot

#if implist, skip imp(), any coordinates, unification, plot
#if missdat, entire work flow


##ggplot

