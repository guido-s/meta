\name{metareg}

\alias{metareg}

\title{Meta-regression}

\description{
  Meta-regression for objects of class \code{meta}. This is a wrapper
  function for the R function \code{\link[metafor]{rma.uni}} in the R package
  metafor (Viechtbauer 2010).
}

\usage{
metareg(x, formula, method.tau=x$method.tau, hakn=x$hakn,
        level.comb=x$level.comb, intercept=TRUE,...)
}

\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{formula}{Either a character string or a formula object.}
  \item{method.tau}{A character string indicating which method is used
    to estimate the between-study variance tau-squared. Either
    \code{"FE"}, \code{"DL"}, \code{"REML"}, \code{"ML"}, \code{"HS"},
    \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be abbreviated.}
  \item{hakn}{A logical indicating whether the method by Hartung and
    Knapp should be used to adjust test statistics and confidence
    intervals.}
  \item{level.comb}{The level used to calculate confidence intervals for
    parameter estimates in the meta-regression model.}
  \item{intercept}{A logical indicating whether an intercept should be
    included in the meta-regression model.}
  \item{\dots}{Additional arguments passed to R function
    \code{\link[metafor]{rma.uni}}.}
}

\details{
  This R function is a wrapper function for R function
  \code{\link[metafor]{rma.uni}} in the R package metafor (Viechtbauer 2010),
  i.e. the function \code{metareg} can only be used if the R package
  'metafor' is installed.

  Argument '\dots' can be used to pass additional arguments to R
  function \code{\link[metafor]{rma.uni}}. For example, argument \code{control}
  to provide a list of control values for the iterative estimation
  algorithm. See help page of R function \code{\link[metafor]{rma.uni}} for
  more details.
}

\value{
An object of class \code{c("metareg", "rma.uni","rma")}. Please look
at the help page of R function \code{rma.uni} for more details on the
output from this function.

In addition, a list \code{.meta} is added to the output containing the
following components:
  \item{x, formula, method.tau, hakn, level.comb, intercept}{As definied above.}
  \item{dots}{Information provided in argument '\dots'.}
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

# Add some (fictious) grouping variables:
Fleiss93cont$age <- c(55, 65, 55, 65, 55)
Fleiss93cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")

meta1 <- metacont(n.e, mean.e, sd.e,
                  n.c, mean.c, sd.c,
                  data=Fleiss93cont, sm="MD")


mu1 <- update(meta1, byvar=region)

mu2 <- update(meta1, byvar=region,
              tau.common=TRUE, comb.fixed=FALSE)

\dontrun{
# Warnings due to wrong ordering of arguments (order has changed with
# version 3.0-0 of R package meta)
#
metareg(~region, meta1)
metareg(~region, data=meta1)

# Warning as no information on covariate is available
#
metareg(meta1)
}

# Do meta-regression for covariate region
# (see R code to create object mu2)
#
metareg(mu2)

# Same result for
# - tau-squared
# - test of heterogeneity
# - test for subgroup differences
# (as argument 'tau.common' was used to create mu2)
#
mu2
metareg(mu2, intercept=FALSE)
metareg(meta1, region)
#
# Different result for
# - tau-squared
# - test of heterogeneity
# - test for subgroup differences
# (as argument 'tau.common' is - by default - FALSE)
#
mu1

# Generate bubble plot
#
bubble(metareg(mu2))

# Do meta-regression with two covariates
#
metareg(mu1, region + age)

# Do same meta-regressions using 'official' formula notation
#
metareg(meta1, ~region)
metareg(mu1, ~region + age)

# Do meta-regression using REML method and print intermediate results
# for iterative estimation algorithm; furthermore print results with
# three digits.
#
metareg(mu1, region, method.tau="REML",
        control=list(verbose=TRUE), digits=3)

# Use Hartung-Knapp method
#
mu3 <- update(mu2, hakn=TRUE)
mu3
metareg(mu3, intercept=FALSE)
}

\keyword{print}
