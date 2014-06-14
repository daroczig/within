# Within estimators for fixed effects panel data

This package includes R functions to the most frequently used three-dimensional fixed effects panel data models and the appropriate Within estimators based on

 * László Balázsi, László Mátyás and Tom Wansbeek (2014): [The Estimation of Multi-dimensional Fixed Effects Panel Data Models](http://EconPapers.repec.org/RePEc:ceu:econwp:2014_1). No 2014_1, CEU Working Papers, Department of Economics, Central European University.

The original goal of this package is to provide a way to do all the required transformations and estimates for reasonable large datasets on commodity hardware in reasonable time. The algorithms were tested on a laptop with 8 Gb of RAM for 180 x 180 x 46 matrices.

## Installation

The package is not on CRAN yet, but he most recent version is hosted on [GitHub](https://github.com/rapporter/rapport), so installation is rather easy with the nifty function from `devtools` package:

```
library(devtools)
install_github('within', daroczig')
```

Or download the [sources in a zip file](https://github.com/daroczig/within/archive/master.zip) and build manually. To do so, please unzip the file to an empty dir and run the following commands there:

```
R CMD build rapport
R CMD INSTALL within_0.1.tar.gz
```

If you're running R on Windows, you need to install [Rtools](http://cran.stat.ucla.edu/bin/windows/Rtools/ ). Once you have installed `Rtools`, issue following command in command prompt:

```
R CMD build --binary <path to .tar.gz file>
R CMD INSTALL <path to .zip file>
```
