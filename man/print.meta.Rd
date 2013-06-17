\name{print.meta}

\alias{print.meta}
\alias{print.metabias}
\alias{summary.meta}
\alias{print.summary.meta}

\title{Print and summary method for objects of class meta}

\description{
  Print and summary method for objects of class \code{meta}.
}

\usage{
\method{print}{meta}(x, sortvar, level=x$level, level.comb=x$level.comb,
        comb.fixed=x$comb.fixed, comb.random=x$comb.random,
        details=FALSE, ma=TRUE, digits=max(4, .Options$digits - 3), ...)

\method{print}{metabias}(x, ...)

\method{summary}{meta}(object, byvar=object$byvar,
        bylab=object$bylab, print.byvar=object$print.byvar,
        bystud=FALSE,
        level=object$level, level.comb=object$level.comb,
        comb.fixed=object$comb.fixed, comb.random=object$comb.random,
        prediction=object$prediction, level.predict=object$level.predict,
        print.CMH=object$print.CMH, warn=object$warn, ...)

\method{print}{summary.meta}(x, digits = max(3, .Options$digits - 3),
        print.byvar=x$print.byvar,
        comb.fixed=x$comb.fixed, comb.random=x$comb.random,
        prediction=x$prediction,
        header=TRUE, print.CMH=x$print.CMH, bylab.nchar=35, ...)
}

\arguments{
  \item{x}{An object of class \code{meta}, \code{metabias}, or
    \code{summary.meta}.}
  \item{object}{An object of class \code{meta}.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{prediction}{A logical indicating whether a prediction interval
    should be printed.}
  \item{level.predict}{The level used to calculate prediction interval
    for a new study.}
  \item{header}{A logical indicating whether information on title of
    meta-analysis, comparison and outcome should be printed at the
    beginning of the printout.}
  \item{details}{A logical indicating whether further details of
    individual studies should be printed.}
  \item{ma}{A logical indicating whether the summary results of the
    meta-analysis should be printed.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{x$TE}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{bylab.nchar}{A numeric specifying the number of characters to
    print from label for the grouping variable.}
  \item{bystud}{A logical indicating whether results of individual
    studies should be printed by grouping variable.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels. By default,
    the value of \code{print.byvar} is set to \code{TRUE}.}
  \item{warn}{A logical indicating whether the use of
    \code{summary.meta} in connection with \code{metacum} or
    \code{metainf} should result in a warning.}
  \item{print.CMH}{A logical indicating whether result of the
    Cochran-Mantel-Haenszel test for overall effect should be printed.}
  \item{...}{Additional arguments}
}

\details{
  Review Manager 5 (RevMan 5) is the current software used for
  preparing and maintaining Cochrane Reviews
  (\url{http://www.cc-ims.net/revman/}). In RevMan 5, subgroup analyses
  can be defined and data from a Cochrane review can be imported to R
  using the function \code{read.rm5}. If a meta-analysis is then
  conducted using function \code{metacr}, information on subgroups is
  available in R (components \code{byvar}, \code{bylab}, and
  \code{print.byvar}, \code{byvar} in an object of class
  \code{"meta"}). Accordingly, by using function \code{metacr} there is
  no need to define subgroups in order to redo the statistical analysis
  conducted in the Cochrane review.

  For subgroups (argument \code{byvar} not \code{NULL}), results for the
  fixed effect model will be printed if both arguments \code{comb.fixed}
  and \code{comb.random} are \code{TRUE}. In order to get results for
  the random effects model within subgroups, use
  \code{comb.fixed==FALSE} and \code{comb.random==TRUE}.

  Note, for an object of type \code{metaprop}, exact binomial
  confidence intervals are calculated for individual study results
  using the R function \code{\link{binom.test}}
  internally. Accordingly, list elements \code{TE}, \code{lower} and
  \code{upper} in element \code{study} correspond to proportions and
  exact confidence limits on the natural scale (irrespective of the
  transformation used in meta-analysis). Contrary, meta-analysis
  results are transformed as defined by argument \code{sm}, i.e. list
  elements \code{TE}, \code{lower} and \code{upper} in elements
  \code{fixed}, \code{random}, \code{within.fixed} and
  \code{within.random}.
}

\value{
  A list is returned by the function \code{summary.meta} with the
  following elements:
  \item{study}{Results for individual studies (a list with elements TE,
    seTE, lower, upper, z, p, level, df).}
  \item{fixed}{Results for fixed effect model (a list with elements TE,
    seTE, lower, upper, z, p, level, df).}
  \item{random}{Results for random effects model (a list with elements
    TE, seTE, lower, upper, z, p, level, df).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance.}
  \item{se.tau}{Standard error of square-root of between-study variance.}
  \item{C}{Scaling factor utilised internally to calculate common
    tau-squared across subgroups.}
  \item{H}{Heterogeneity statistic H  (a list with elements TE, lower, upper).}
  \item{I2}{Heterogeneity statistic I2 (a list with elements TE, lower,
    upper), see Higgins & Thompson (2002).}
  \item{k.all}{Total number of studies.}
  \item{Q.CMH}{Cochran-Mantel-Haenszel test statistic for overall effect.}
  \item{sm}{A character string indicating underlying summary measure.}
  \item{method}{A character string with the pooling method.}
  \item{call}{Function call.}
  \item{ci.lab}{Label for confidence interval.}
  \item{hakn}{A logical indicating whether method by Hartung and
    Knapp was used.}
  \item{method.tau}{A character string indicating which method is used
    to estimate the between-study variance tau-squared.}
  \item{tau.common}{A logical indicating whether tau-squared is
    assumed to be the same across subgroups.}
  \item{within.fixed}{Result for fixed effect model within groups (a
    list with elements TE, seTE, lower, upper, z, p, level, df,
    harmonic.mean) - if \code{byvar} is not missing.}
  \item{within.random}{Result for random effects model within groups
    (a list with elements TE, seTE, lower, upper, z, p, level, df,
    harmonic.mean) - if \code{byvar} is not missing.}
  \item{k.w}{Number of studies combined within groups - if \code{byvar}
    is not missing.}
  \item{Q.w}{Heterogeneity statistic Q within groups - if \code{byvar}
    is not missing.}
  \item{Q.b.fixed}{Heterogeneity statistic Q between groups (based on
    fixed effect model) - if \code{byvar} is not missing.}
  \item{Q.b.random}{Heterogeneity statistic Q between groups (based on
    random effects model) - if \code{byvar} is not missing.}
  \item{tau.w}{Square-root of between-study variance within subgroups
    - if \code{byvar} is not missing.}
  \item{C.w}{Scaling factor utilised internally to calculate common
    tau-squared across subgroups.}
  \item{H.w}{Heterogeneity statistic H within subgroups (a list with
    elements TE, lower, upper) - if \code{byvar} is not missing.}
  \item{I2.w}{Heterogeneity statistic I2 within subgroups (a list with
    elements TE, lower, upper) - if \code{byvar} is not missing.}
  \item{bylab}{Label for grouping variable - if \code{byvar} is not
    missing.}
  \item{by.levs}{Levels of grouping variable - if \code{byvar} is not
    missing.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{prediction, level.predict}{}
  \item{comb.fixed, comb.random, print.CMH}{As defined above.}
  \item{version}{Version of R package meta used to create object.}
}

\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.

  Higgins JPT & Thompson SG (2002),
  Quantifying heterogeneity in a meta-analysis.
  \emph{Statistics in Medicine}, \bold{21}, 1539--1558.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
data(Fleiss93cont)
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont, sm="SMD")
summary(meta1)
summary(meta1, byvar=c(1,2,1,1,2), bylab="group")
}

\keyword{print}
