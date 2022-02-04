# meta: General Package for Meta-Analysis
Official Git repository of R package **meta**

[![Build Status](https://travis-ci.org/guido-s/meta.svg?branch=master)](https://travis-ci.org/guido-s/meta)
[![CRAN Version](http://www.r-pkg.org/badges/version/meta)](https://cran.r-project.org/package=meta)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/meta)](http://cranlogs.r-pkg.org/badges/meta)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/meta)](http://cranlogs.r-pkg.org/badges/grand-total/meta)


## Description

User-friendly general package providing standard methods for meta-analysis and supporting Schwarzer, Carpenter, and Rücker, "Meta-Analysis with R" (2015):
 - fixed effect and random effects meta-analysis;
 - several plots (forest, funnel, Galbraith / radial, L'Abbe, Baujat, bubble);
 - statistical tests and trim-and-fill method to evaluate bias in meta-analysis;
 - import data from 'RevMan 5';
 - prediction interval, Hartung-Knapp method for random effects model;
 - cumulative meta-analysis and leave-one-out meta-analysis;
 - meta-regression;
 - generalised linear mixed models;
 - produce forest plot summarising several (subgroup) meta-analyses.
 
### References

[Schwarzer G, Carpenter JR and Rücker G (2015): *Meta-Analysis with R (Use-R!)*. Springer International Publishing, Switzerland](https://link.springer.com/book/10.1007/978-3-319-21416-0)


## Installation

### Current official [![CRAN Version](http://www.r-pkg.org/badges/version/meta)](https://cran.r-project.org/package=meta) release:
```r
install.packages("meta")
```

### Current beta / GitHub release:

Installation using R package
[**remotes**](https://cran.r-project.org/package=remotes):
```r
install.packages("remotes")
remotes::install_github("guido-s/meta", ref = "develop")
```


### Bug Reports:

```r
bug.report(package = "meta")
```

The bug.report function is not supported in RStudio. Please send an
email to Guido Schwarzer <sc@imbi.uni-freiburg.de> if you use RStudio.

You can also report bugs on GitHub under [Issues](https://github.com/guido-s/meta/issues/).
