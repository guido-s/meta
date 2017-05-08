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
        pscale=x$pscale, irscale = x$irscale, irunit = x$irunit,
        digits = gs("digits"), digits.se = gs("digits.se"),
        digits.zval = gs("digits.zval"),
        digits.pval = max(gs("digits.pval"), 2),
        digits.pval.Q = max(gs("digits.pval.Q"), 2),
        digits.Q = gs("digits.Q"), digits.tau2 = gs("digits.tau2"),
        digits.H = gs("digits.H"), digits.I2 = gs("digits.I2"),
        digits.prop = gs("digits.prop"), digits.weight = gs("digits.weight"),
	scientific.pval = gs("scientific.pval"),
        warn.backtransf = FALSE,
        ...)

\method{print}{metabias}(x, ...)

\method{summary}{meta}(object,
        comb.fixed=object$comb.fixed, comb.random=object$comb.random,
        prediction=object$prediction,
        backtransf=object$backtransf,
        pscale=object$pscale, irscale = object$irscale, irunit = object$irunit,
        bylab=object$bylab, print.byvar=object$print.byvar,
        byseparator=object$byseparator, bystud=FALSE,
        print.CMH=object$print.CMH, warn=object$warn, ...)

\method{print}{summary.meta}(x, digits = gs("digits"),
        comb.fixed=x$comb.fixed, comb.random=x$comb.random,
        prediction=x$prediction,
        print.byvar=x$print.byvar, byseparator=x$byseparator,
        print.CMH=x$print.CMH,
        header=TRUE, backtransf=x$backtransf,
        pscale=x$pscale, irscale = x$irscale, irunit = x$irunit,
        bylab.nchar=35,
        digits.zval = gs("digits.zval"),
        digits.pval = max(gs("digits.pval"), 2),
        digits.pval.Q = max(gs("digits.pval.Q"), 2),
        digits.Q = gs("digits.Q"), digits.tau2 = gs("digits.tau2"),
        digits.H = gs("digits.H"), digits.I2 = gs("digits.I2"),
	scientific.pval = gs("scientific.pval"),
        print.I2 = gs("print.I2"), print.H = gs("print.H"),
        print.Rb = gs("print.Rb"),
        text.tau2 = gs("text.tau2"), text.I2 = gs("text.I2"),
        text.Rb = gs("text.Rb"),
        warn.backtransf = FALSE,
        ...)

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
  \item{byseparator}{A character string defining the separator between
    label and levels of grouping variable.}
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
  \item{pscale}{A numeric giving scaling factor for printing of single
    event probabilities, i.e. if argument \code{sm} is equal to
    \code{"PLOGIT"}, \code{"PLN"}, \code{"PRAW"}, \code{"PAS"}, or
    \code{"PFT"}.}
  \item{irscale}{A numeric defining a scaling factor for printing of
    rates, i.e. if argument \code{sm} is equal to \code{"IR"},
    \code{"IRLN"}, \code{"IRS"}, or \code{"IRFT"}.}
  \item{irunit}{A character specifying the time unit used to calculate
    rates, e.g. person-years.}
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
  \item{warn.backtransf}{A logical indicating whether a warning should
    be printed if backtransformed proportions and rates are below 0
    and backtransformed proportions are above 1.}
  \item{bracket}{A character with bracket symbol to print lower
    confidence interval: "[", "(", "\{", "".}
  \item{separator}{A character string with information on separator
    between lower and upper confidence interval.}
  \item{digits.se}{Minimal number of significant digits for standard
    deviations and standard errors, see \code{print.default}.}
  \item{digits.zval}{Minimal number of significant digits for z- or
    t-value, see \code{print.default}.}
  \item{digits.pval}{Minimal number of significant digits for p-value
    of overall treatment effect, see \code{print.default}.}
  \item{digits.pval.Q}{Minimal number of significant digits for
    p-value of heterogeneity test, see \code{print.default}.}
  \item{digits.Q}{Minimal number of significant digits for
    heterogeneity statistic Q, see \code{print.default}.}
  \item{digits.tau2}{Minimal number of significant digits for
    between-study variance, see \code{print.default}.}
  \item{digits.H}{Minimal number of significant digits for H
    statistic, see \code{print.default}.}
  \item{digits.I2}{Minimal number of significant digits for I-squared
    and Rb statistic, see \code{print.default}.}
  \item{digits.prop}{Minimal number of significant digits for
    proportions, see \code{print.default}.} 
  \item{digits.weight}{Minimal number of significant digits for
    weights, see \code{print.default}.}    
  \item{scientific.pval}{A logical specifying whether p-values should
    be printed in scientific notation, e.g., 1.2345e-01 instead of
    0.12345.}
  \item{print.I2}{A logical specifying whether heterogeneity statistic
    I^2 should be printed.}
  \item{print.H}{A logical specifying whether heterogeneity statistic
    H should be printed.}
  \item{print.Rb}{A logical specifying whether heterogeneity statistic
    Rb should be printed.}
  \item{text.tau2}{Text printed to identify between-study variance tau^2.}
  \item{text.I2}{Text printed to identify heterogeneity statistic I^2.}
  \item{text.Rb}{Text printed to identify heterogeneity statistic Rb.}
  \item{\dots}{Additional arguments.}
}

\details{
  Note, in R package \bold{meta}, version 3.0-0 some arguments have
  been removed from R functions \code{\link{summary.meta}} (arguments:
  byvar, level, level.comb, level.prediction) and print.summary.meta
  (arguments: level, level.comb, level.prediction). This functionality
  is now provided by R function \code{\link{update.meta}} (or directly
  in meta-analysis functions, e.g., \code{\link{metabin}},
  \code{\link{metacont}}, \code{\link{metagen}},
  \code{\link{metacor}}, and \code{\link{metaprop}}).

  Review Manager 5 (RevMan 5) is the current software used for
  preparing and maintaining Cochrane Reviews
  (\url{http://community.cochrane.org/tools/review-production-tools/revman-5}). In
  RevMan 5, subgroup analyses can be defined and data from a Cochrane
  review can be imported to R using the function \code{read.rm5}. If a
  meta-analysis is then conducted using function \code{metacr},
  information on subgroups is available in R (components \code{byvar},
  \code{bylab}, and \code{print.byvar}, \code{byvar} in an object of
  class \code{"meta"}). Accordingly, by using function \code{metacr}
  there is no need to define subgroups in order to redo the
  statistical analysis conducted in the Cochrane review.
  
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
  
  Argument \code{pscale} can be used to rescale proportions,
  e.g. \code{pscale=1000} means that proportions are expressed as
  events per 1000 observations. This is useful in situations with
  (very) low event probabilities.
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
  \item{Rb}{Heterogeneity statistic Rb (a list with elements TE,
    lower, upper), see Crippa et al. (2016).}
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
  \item{Rb.w}{Heterogeneity statistic Rb within subgroups (a list with
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
  \item{version}{Version of R package \bold{meta} used to create object.}
}

\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.

  Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D (2016),
  A new measure of between-studies heterogeneity in meta-analysis.
  \emph{Statistics in Medicine}, \bold{35}, 3661--75.

  Higgins JPT & Thompson SG (2002),
  Quantifying heterogeneity in a meta-analysis.
  \emph{Statistics in Medicine}, \bold{21}, 1539--58.
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

\dontrun{
# Use unicode characters to print tau^2 and I^2 
#
print(summary(meta1), text.tau2 = "\u03c4\u00b2", text.I2 = "I\u00b2")
}
}

\keyword{print}
