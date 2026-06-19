# meta: General Package for Meta-Analysis
Official Git repository of R package **meta**

[![License: GPL (>=2)](https://img.shields.io/badge/license-GPL-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![CRAN Version](https://www.r-pkg.org/badges/version/meta)](https://cran.r-project.org/package=meta)
[![GitHub develop](https://img.shields.io/badge/develop-8.5--1-purple)](https://img.shields.io/badge/develop-8.5--1-purple)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/meta)](https://cranlogs.r-pkg.org/badges/meta)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/meta)](https://cranlogs.r-pkg.org/badges/grand-total/meta)


## Author

[Guido Schwarzer](https://orcid.org/0000-0001-6214-9087)


## Description

User-friendly general package providing standard methods for meta-analysis and supporting Schwarzer, Carpenter, and Rücker, "Meta-Analysis with R" (2015):
 - common effect and random effects meta-analysis;
 - several plots (forest, funnel, Galbraith / radial, L'Abbe, Baujat, bubble);
 - statistical tests and trim-and-fill method to evaluate bias in meta-analysis;
 - Hartung-Knapp method for random effects model;
 - prediction interval and density of the prediction distribution;
 - expected proportion of comparable studies with clinically important benefit or harm;
 - cumulative meta-analysis and leave-one-out meta-analysis;
 - meta-regression;
 - generalised linear mixed models;
 - logistic regression with penalised likelihood for rare events;
 - import data from 'RevMan Web' and 'RevMan 5';
 - produce forest plot summarising several (subgroup) meta-analyses.
 
### References

[Schwarzer G, Carpenter JR and Rücker G (2015): *Meta-Analysis with R (Use-R!)*. Springer International Publishing, Switzerland](https://link.springer.com/book/10.1007/978-3-319-21416-0)


## Installation

### Current official [![CRAN Version](https://www.r-pkg.org/badges/version/meta)](https://cran.r-project.org/package=meta) release:
```r
install.packages("meta")
```

### Current [![GitHub develop](https://img.shields.io/badge/develop-8.5--1-purple)](https://img.shields.io/badge/develop-8.5--1-purple) release on GitHub:

Installation using R package
[**remotes**](https://cran.r-project.org/package=remotes):
```r
install.packages("remotes")
remotes::install_github("guido-s/meta", build_vignettes = TRUE)
```
or without the vignette
```r
remotes::install_github("guido-s/meta")
```


## How to cite netmeta?

[Balduzzi S, Rücker G, Schwarzer G (2019): How to perform a meta-analysis with R: a practical tutorial. *Evidence-Based Mental Health*, **22**, 153-60](https://doi.org/10.1136/ebmental-2019-300117)

A BibTeX entry for LaTeX users is provided by

```
citation(package = "meta")
```


### Bug Reports:

You can report bugs on GitHub under
[Issues](https://github.com/guido-s/meta/issues)

or by using the R command

```r
bug.report(package = "meta")
```

(which is not supported in RStudio).
