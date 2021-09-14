# meta: General Package for Meta-Analysis
Git repository of R package **meta** to support [*Meta-Analysis with R
(Use-R!), first
edition*](http://www.springer.com/gp/book/9783319214153)


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

[Schwarzer G, Carpenter JR and Rücker G (2015): *Meta-Analysis with R (Use-R!)*. Springer International Publishing, Switzerland](http://www.springer.com/gp/book/9783319214153)


## Installation

Installation using R package
[**devtools**](https://cran.r-project.org/package=devtools):
```r
install.packages("devtools")
devtools::install_github("guido-s/meta, ref = "R-book-first-edition")
```

Please send an email to Guido Schwarzer <sc@imbi.uni-freiburg.de> for
additional support on this version of **meta**.
