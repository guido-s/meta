\name{metareg}

\alias{metareg}

\title{Meta-regression}

\description{
  Meta-regression for objects of class \code{meta}. This is a wrapper
  function for the R function \code{rma.uni} in the R package metafor
  (Viechtbauer 2010).
}

\usage{
metareg(x, formula, method.tau=x$method.tau, ...)
}

\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{formula}{Either a character string or a formula object.}
  \item{method.tau}{A character string indicating which method is used
    to estimate the between-study variance tau-squared. Either
    \code{"FE"}, \code{"DL"}, \code{"REML"}, \code{"ML"}, \code{"HS"},
    \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be abbreviated.}
  \item{\dots}{Additional arguments (ignored at the moment).}
}

\details{
  This R function is a wrapper function for R function \code{rma.uni}
  in the R package metafor (Viechtbauer 2010), i.e. the function
  \code{metareg} can only be used if the R package metafor is
  installed.
}

\value{
An object of class \code{c("metareg", "rma.uni","rma")}. Please look
at the help page of R function \code{rma.uni} for more details on the
output from this function.

In addition, a list \code{.meta} is added to the output containing the
following components:
  \item{x, formula, method.tau}{As definied above.}
  \item{call}{Function call.}
  \item{version}{Version of R package meta used to create object.}
  \item{version.metafor}{Version of R package metafor used to create object.}
}

\references{
  Viechtbauer W (2010),
  Conducting Meta-Analyses in R with the Metafor Package.
  \emph{Journal of Statistical Software}, \bold{36}, 1--48.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{bubble}}, \code{\link{summary.meta}}, \code{\link{metagen}}}

\examples{
data(Fleiss93cont)

## Add some (fictious) grouping variables:
Fleiss93cont$age <- c(55, 65, 55, 65, 55)
Fleiss93cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")

meta1 <- metacont(n.e, mean.e, sd.e,
                  n.c, mean.c, sd.c,
                  data=Fleiss93cont, sm="MD")


mu1 <- update(meta1, byvar=region)

mu2 <- update(meta1, byvar=region,
              tau.common=TRUE, comb.fixed=FALSE)

\dontrun{
## Warnings due to wrong ordering of arguments (order has changed with
## version 3.0-0 of R package meta)
##
metareg(~region, meta1)
metareg(~region, data=meta1)

## Warning as no information on covariate is available
##
metareg(meta1)
}

## Do meta-regression for covariate region
## (see R code to create object mu2)
##
metareg(mu2)

## Same result for
## - tau-squared
## - test of heterogeneity
## - test for subgroup differences
## (as argument 'tau.common' was used to create mu2)
##
mu2
metareg(mu2)
metareg(meta1, region)
##
## Different result for
## - tau-squared
## - test of heterogeneity
## - test for subgroup differences
## (as argument 'tau.common' is - by default - FALSE)
##
mu1

## Generate bubble plot
##
bubble(metareg(mu2))

## Do meta-regression with two covariates
##
metareg(mu1, region + age)

## Do same meta-regressions using 'official' formula notation
##
metareg(meta1, ~region)
metareg(mu1, ~region + age)
}

\keyword{print}
