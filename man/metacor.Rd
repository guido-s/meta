\name{metacor}
\alias{metacor}

\title{Meta-analysis of correlations}

\description{
  Calculation of fixed and random effects estimates for meta-analyses
  with correlations; inverse variance weighting is used for pooling.
}

\usage{
metacor(cor, n, studlab,
        data=NULL, subset=NULL,
        sm="ZCOR",
        level=0.95, level.comb=level,
        comb.fixed=TRUE, comb.random=TRUE,
        hakn=FALSE,
        method.tau="DL", tau.preset=NULL, TE.tau=NULL,
        tau.common=FALSE,
        prediction=FALSE, level.predict=level,
        method.bias="linreg",
        title="", complab="", outclab="",
        byvar, bylab, print.byvar=TRUE,
        keepdata=TRUE
        )
}

\arguments{
  \item{cor}{Correlation.}
  \item{n}{Number of observations.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., cor and n.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating which summary measure
    (\code{"ZCOR"} or \code{"COR"}) is to be used for pooling of
    studies.}
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
  \item{hakn}{A logical indicating whether the method by Hartung and
          Knapp should be used to adjust test statistics and
          confidence intervals.}
  \item{method.tau}{A character string indicating which method is used
    to estimate the between-study variance tau-squared. Either
    \code{"DL"}, \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
    \code{"HE"}, or \code{"EB"}, can be abbreviated.}
  \item{tau.preset}{Prespecified value for between-study variance tau-squared.}
  \item{TE.tau}{Overall treatment effect used to estimate the
    between-study variance tau-squared.}
  \item{tau.common}{A logical indicating whether tau-squared should be
    the same across subgroups.}
  \item{method.bias}{A character string indicating which test is to be
    used. Either \code{"rank"}, \code{"linreg"}, or \code{"mm"}, can
    be abbreviated.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{event.e}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in meta object.}
}

\details{
  Fixed effect and random effects meta-analysis of correlations
  based either on Fisher's z transformation of correlations (\code{sm="ZCOR"})
  or direct combination of correlations (\code{sm="COR"}) (see Cooper et
  al., p264-5 and p273-4).

  Only few statisticians would advocate the use of untransformed
  correlations unless sample sizes are very large (see Cooper et al.,
  p265). The artificial example given below shows that the smallest
  study gets the largest weight if correlations are combined directly
  because the correlation is closest to 1.

  Internally, both fixed effect and random effects models are calculated
  regardless of values choosen for arguments \code{comb.fixed} and
  \code{comb.random}. Accordingly, the estimate for the random effects
  model can be extracted from component \code{TE.random} of an object
  of class \code{"meta"} even if \code{comb.random=FALSE}. However, all
  functions in R package \code{meta} will adequately consider the values
  for \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.meta}} will not print results for the random effects
  model if \code{comb.random=FALSE}.

  If R package metafor (Viechtbauer 2010) is installed, the following
  statistical methods are also available.

  For the random effects model (argument \code{comb.random=TRUE}), the
  method by Hartung and Knapp (Knapp, Hartung 2003) is used to adjust
  test statistics and confidence intervals if argument
  \code{hakn=TRUE} (internally R function \code{rma.uni} of R package
  metafor is called).

  Several methods are available to estimate the between-study variance
  \eqn{\tau^2} (argument \code{method.tau}):
  \itemize{
   \item DerSimonian-Laird estimator (\code{method.tau="DL"}) (default)
   \item Restricted maximum-likelihood estimator (\code{method.tau="REML"})
   \item Maximum-likelihood estimator (\code{method.tau="ML"})
   \item Hunter-Schmidt estimator (\code{method.tau="HS"})
   \item Sidik-Jonkman estimator (\code{method.tau="SJ"})
   \item Hedges estimator (\code{method.tau="HE"})
   \item Empirical Bayes estimator (\code{method.tau="EB"}).
  }
  For all but the DerSimonian-Laird method the R function
  \code{rma.uni} of R package metafor is called internally. See help
  page of R function \code{rma.uni} for more details on the various
  methods to estimate between-study variance \eqn{\tau^2}{tau^2}.

  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metacor object by only specifying
  arguments which should be changed.
}

\value{
  An object of class \code{c("metacor", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{cor, n, studlab,}{}
  \item{sm, level, level.comb,}{}
  \item{comb.fixed, comb.random,}{}
  \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{}
  \item{tau.common, title, complab, outclab,}{}
  \item{byvar, bylab, print.byvar}{As defined above.}
  \item{TE, seTE}{Either Fisher's z transformation of correlations
    (\code{sm="ZCOR"}) or correlations (\code{sm="COR"}) for individual
    studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall effect (Fisher's z
    transformation of correlation or correlation) and its standard error
    (fixed effect model).}
  \item{lower.fixed, upper.fixed}{Lower and upper confidence interval limits
  (fixed effect model).}
  \item{zval.fixed, pval.fixed}{z-value and p-value for test of
  overall effect (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall effect (Fisher's z
    transformation of correlation or correlation) and its standard error
    (random effects model).}
  \item{lower.random, upper.random}{Lower and upper confidence interval limits
  (random effects model).}
  \item{zval.random, pval.random}{z-value or t-value and corresponding
    p-value for test of overall effect (random effects model).}
  \item{prediction, level.predict}{As defined above.}
  \item{seTE.predict}{Standard error utilised for prediction interval.}
  \item{lower.predict, upper.predict}{Lower and upper limits of prediction interval.}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance.}
  \item{se.tau}{Standard error of square-root of between-study variance.}
  \item{C}{Scaling factor utilised internally to calculate common
    tau-squared across subgroups.}
  \item{method}{A character string indicating method used
    for pooling: \code{"Inverse"}}
  \item{df.hakn}{Degrees of freedom for test of treatment effect for
  Hartung-Knapp method (only if \code{hakn=TRUE}).}
  \item{keepdata}{As defined above.}
  \item{data}{Original data (set) used in function call (if
    \code{keepdata=TRUE}).}
  \item{subset}{Information on subset of original data used in
    meta-analysis (if \code{keepdata=TRUE}).}
  \item{call}{Function call.}
  \item{version}{Version of R package meta used to create object.}
}

\references{
  Cooper H, Hedges LV, Valentine JC (2009),
  \emph{The Handbook of Research Synthesis and Meta-Analysis}, 2nd Edition.
  New York: Russell Sage Foundation.
  
  Higgins JPT, Thompson SG, Spiegelhalter DJ (2009),
  A re-evaluation of random-effects meta-analysis.
  \emph{Journal of the Royal Statistical Society: Series A},
  \bold{172}, 137-159.
  
  Knapp G & Hartung J (2003),
  Improved Tests for a Random Effects Meta-regression with a Single
  Covariate.
  \emph{Statistics in Medicine}, \bold{22}, 2693-710,
  doi: 10.1002/sim.1482 .

  Viechtbauer W (2010),
  Conducting Meta-Analyses in R with the Metafor Package.
  \emph{Journal of Statistical Software}, \bold{36}, 1--48.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{update.meta}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))
forest(metacor(c(0.85, 0.7, 0.95), c(20, 40, 10)))

metacor(c(0.85, 0.7, 0.95), c(20, 40, 10), sm="cor")
}

\keyword{htest}
