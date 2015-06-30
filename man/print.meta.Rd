\name{print.meta}

\alias{print.meta}
\alias{print.metabias}
\alias{summary.meta}
\alias{print.summary.meta}
\alias{cilayout}

\title{Print and summary method for objects of class meta}

\description{
  Print and summary method for objects of class \code{meta}.
}

\usage{
\method{print}{meta}(x, sortvar,
        comb.fixed=x$comb.fixed,
        comb.random=x$comb.random,
        prediction=x$prediction,
        details=FALSE, ma=TRUE, backtransf=x$backtransf,
        digits=max(4, .Options$digits - 3), ...)

\method{print}{metabias}(x, ...)

\method{summary}{meta}(object,
        comb.fixed=object$comb.fixed, comb.random=object$comb.random,
        prediction=object$prediction,
        backtransf=object$backtransf,
        bylab=object$bylab, print.byvar=object$print.byvar,
        bystud=FALSE, print.CMH=object$print.CMH,
        warn=object$warn, ...)

\method{print}{summary.meta}(x, digits = max(3, .Options$digits - 3),
        comb.fixed=x$comb.fixed, comb.random=x$comb.random,
        prediction=x$prediction,
        print.byvar=x$print.byvar, print.CMH=x$print.CMH,
        header=TRUE, backtransf=x$backtransf,
        bylab.nchar=35, ...)

cilayout(bracket="[", separator="; ")
}

\arguments{
  \item{x}{An object of class \code{meta}, \code{metabias}, or
    \code{summary.meta}.}
  \item{object}{An object of class \code{meta}.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{prediction}{A logical indicating whether a prediction interval
    should be printed.}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{header}{A logical indicating whether information on title of
    meta-analysis, comparison and outcome should be printed at the
    beginning of the printout.}
  \item{details}{A logical indicating whether further details of
    individual studies should be printed.}
  \item{ma}{A logical indicating whether the summary results of the
    meta-analysis should be printed.}
  \item{backtransf}{A logical indicating whether printed results
    should be back transformed. If \code{backtransf=TRUE}, results for
    \code{sm="OR"} are printed as odds ratios rather than log odds
    ratios and results for \code{sm="ZCOR"} are printed as
    correlations rather than Fisher's z transformed correlations, for
    example.}
  \item{bylab.nchar}{A numeric specifying the number of characters to
    print from label for the grouping variable.}
  \item{bystud}{A logical indicating whether results of individual
    studies should be printed by grouping variable.}
  \item{print.CMH}{A logical indicating whether result of the
    Cochran-Mantel-Haenszel test for overall effect should be printed.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{warn}{A logical indicating whether the use of
    \code{summary.meta} in connection with \code{metacum} or
    \code{metainf} should result in a warning.}
  \item{bracket}{A character with bracket symbol to print lower
    confidence interval: "[", "(", "\{", "".}
  \item{separator}{A character string with information on separator
    between lower and upper confidence interval.}
  \item{\dots}{Additional arguments.}
}

\details{
  Note, in R package meta, version 3.0-0 some arguments have been
  removed from R functions summary.meta (arguments: byvar, level,
  level.comb, level.prediction) and print.summary.meta (arguments:
  level, level.comb, level.prediction). This functionality is now
  provided by R function \code{\link{update.meta}} (or directly in R
  functions \code{\link{metabin}}, \code{\link{metacont}},
  \code{\link{metagen}}, \code{\link{metacor}}, and
  \code{\link{metaprop}}).

  Review Manager 5 (RevMan 5) is the current software used for
  preparing and maintaining Cochrane Reviews
  (\url{http://tech.cochrane.org/revman/}). In RevMan 5, subgroup analyses
  can be defined and data from a Cochrane review can be imported to R
  using the function \code{read.rm5}. If a meta-analysis is then
  conducted using function \code{metacr}, information on subgroups is
  available in R (components \code{byvar}, \code{bylab}, and
  \code{print.byvar}, \code{byvar} in an object of class
  \code{"meta"}). Accordingly, by using function \code{metacr} there is
  no need to define subgroups in order to redo the statistical analysis
  conducted in the Cochrane review.
  
  Note, for an object of type \code{metaprop}, starting with version
  3.7-0 of meta, list elements \code{TE}, \code{lower} and
  \code{upper} in element \code{study} correspond to transformed
  proportions and confidence limits (regardless whether exact
  confidence limits are calculated; argument \code{ciexact=TRUE} in
  metaprop function). Accordingly, the following results are based on
  the same transformation defined by argument \code{sm}: list elements
  \code{TE}, \code{lower} and \code{upper} in elements \code{study},
  \code{fixed}, \code{random}, \code{within.fixed} and
  \code{within.random}.

  R function cilayout can be utilised to change the layout to print
  confidence intervals (both in printout from print.meta and
  print.summary.meta function as well as in forest plots). The default
  layout is "[lower; upper]". Another popular layout is "(lower -
  upper)" which is used throughout an R session by using R command
  \code{cilayout("(", " - ")}.
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
  \item{bylevs}{Levels of grouping variable - if \code{byvar} is not
    missing.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{data}{Original data (set) used to create meta object.}
  \item{subset}{Information on subset of original data used in meta-analysis.}
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

\seealso{\code{\link{update.meta}}, \code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
data(Fleiss93cont)
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
                  data=Fleiss93cont, sm="SMD",
                  studlab=paste(study, year))
summary(meta1)

summary(update(meta1, byvar=c(1,2,1,1,2), bylab="group"))

forest(update(meta1, byvar=c(1,2,1,1,2), bylab="group"))
}

\keyword{print}
