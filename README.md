# GPAbin
GPAbin biplots to align multidimensional visualisations from multiple imputations.

In order to effortlessly apply the presented methodology in "GPAbin: unifying visualizations of multiple mputations for missing values" in **Communications in Statistics: Simulation and Computation**. 52:6, 2666-2685. (https://doi.org/10.1080/03610918.2021.1914089)

Authors: J Nienkemper-Swanepoel, NJ le Roux & S Gardner-Lubbe
Centre for Multi-dimensional data visualisation (MuViSU), Department of Statistics and Actuarial Science, Stellenbosch University.

To start:
Open Code > GPAbincalls.R
Run the function calls in line 30-63.

Users have the option of utilising an example data set, *comp.dat.txt*, which is a fully observed simulated data set (uniform distribution n=100, p=5). This is accompanied by *miss.dat.txt* which is the incomplete version of *comp.dat.txt* with 10% missing values inserted with a missing at random mechanism.

Alternatively, incomplete categorical multivariate data sets may be inserted in the place of object *miss.dat*.
