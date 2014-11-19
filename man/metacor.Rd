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
        sm=.settings$smcor,
        level=.settings$level, level.comb=.settings$level.comb,
        comb.fixed=.settings$comb.fixed, comb.random=.settings$comb.random,
        hakn=.settings$hakn,
        method.tau=.settings$method.tau, tau.preset=NULL, TE.tau=NULL,
        tau.common=.settings$tau.common,
        prediction=.settings$prediction, level.predict=.settings$level.predict,
        method.bias=.settings$method.bias,
        backtransf=.settings$backtransf,
        title=.settings$title, complab=.settings$complab, outclab="",
        byvar, bylab, print.byvar=.settings$print.byvar,
        keepdata=.settings$keepdata
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
    to estimate the between-study variance \eqn{\tau^2}. Either
    \code{"DL"}, \code{"PM"}, \code{"REML"}, \code{"ML"}, \code{"HS"},
    \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be abbreviated.}
  \item{tau.preset}{Prespecified value for the square-root of the
    between-study variance \eqn{\tau^2}.}
  \item{TE.tau}{Overall treatment effect used to estimate the
    between-study variance tau-squared.}
  \item{tau.common}{A logical indicating whether tau-squared should be
    the same across subgroups.}
  \item{method.bias}{A character string indicating which test is to be
    used. Either \code{"rank"}, \code{"linreg"}, or \code{"mm"}, can
    be abbreviated. See function \code{\link{metabias}}}
  \item{backtransf}{A logical indicating whether results for Fisher's
    z transformed correlations (\code{sm="ZCOR"}) should be back
    transformed in printouts and plots. If TRUE (default), results
    will be presented as correlations; otherwise Fisher's z
    transformed correlations will be shown.}
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
  Fixed effect and random effects meta-analysis of correlations based
  either on Fisher's z transformation of correlations
  (\code{sm="ZCOR"}) or direct combination of correlations
  (\code{sm="COR"}) (see Cooper et al., p264-5 and p273-4). By
  default, the DerSimonian-Laird estimate (1986) is used in the random
  effects model (\code{method.tau="DL"}).

  Only few statisticians would advocate the use of untransformed
  correlations unless sample sizes are very large (see Cooper et al.,
  p265). The artificial example given below shows that the smallest
  study gets the largest weight if correlations are combined directly
  because the correlation is closest to 1.
  
  For several arguments defaults settings are utilised (assignments
  with \code{.settings$}). These defaults can be changed using the
  \code{\link{settings.meta}} function.
  
  Internally, both fixed effect and random effects models are
  calculated regardless of values choosen for arguments
  \code{comb.fixed} and \code{comb.random}. Accordingly, the estimate
  for the random effects model can be extracted from component
  \code{TE.random} of an object of class \code{"meta"} even if
  argument \code{comb.random=FALSE}. However, all functions in R
  package \code{meta} will adequately consider the values for
  \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.meta}} will not print results for the random
  effects model if \code{comb.random=FALSE}.

  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metacor object by only specifying
  arguments which should be changed.

  For the random effects, the method by Hartung and Knapp (2003) is
  used to adjust test statistics and confidence intervals if argument
  \code{hakn=TRUE}.
  
  The iterative Paule-Mandel method (1982) to estimate the
  between-study variance is used if argument
  \code{method.tau="PM"}. Internally, R function \code{paulemandel} is
  called which is based on R function mpaule.default from R package
  metRology from S.L.R. Ellison <s.ellison at lgc.co.uk>.

  If R package metafor (Viechtbauer 2010) is installed, the following
  methods to estimate the between-study variance \eqn{\tau^2}
  (argument \code{method.tau}) are also available:
  \itemize{
   \item Restricted maximum-likelihood estimator (\code{method.tau="REML"})
   \item Maximum-likelihood estimator (\code{method.tau="ML"})
   \item Hunter-Schmidt estimator (\code{method.tau="HS"})
   \item Sidik-Jonkman estimator (\code{method.tau="SJ"})
   \item Hedges estimator (\code{method.tau="HE"})
   \item Empirical Bayes estimator (\code{method.tau="EB"}).
  }
  For these methods the R function \code{rma.uni} of R package metafor
  is called internally. See help page of R function \code{rma.uni} for
  more details on these methods to estimate between-study variance.
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
  \item{lower, upper}{Lower and upper confidence interval limits
    for individual studies.}
  \item{zval, pval}{z-value and p-value for test of treatment
    effect for individual studies.}
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
  \item{bylevs}{Levels of grouping variable - if \code{byvar} is not
    missing.}
  \item{TE.fixed.w, seTE.fixed.w}{Estimated treatment effect and
    standard error in subgroups (fixed effect model) - if \code{byvar}
    is not missing.}
  \item{lower.fixed.w, upper.fixed.w}{Lower and upper confidence
    interval limits in subgroups (fixed effect model) - if
    \code{byvar} is not missing.}
  \item{zval.fixed.w, pval.fixed.w}{z-value and p-value for test of
    treatment effect in subgroups (fixed effect model) - if
    \code{byvar} is not missing.}
  \item{TE.random.w, seTE.random.w}{Estimated treatment effect and
    standard error in subgroups (random effects model) - if
    \code{byvar} is not missing.}
  \item{lower.random.w, upper.random.w}{Lower and upper confidence
    interval limits in subgroups (random effects model) - if
    \code{byvar} is not missing.}
  \item{zval.random.w, pval.random.w}{z-value or t-value and
    corresponding p-value for test of treatment effect in subgroups
    (random effects model) - if \code{byvar} is not missing.}
  \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
    random effects model) - if \code{byvar} is not missing.}
  \item{df.hakn.w}{Degrees of freedom for test of treatment effect for
    Hartung-Knapp method in subgroups - if \code{byvar} is not missing
    and \code{hakn=TRUE}.}
  \item{n.harmonic.mean.w}{Harmonic mean of number of observations in
    subgroups (for back transformation of Freeman-Tukey Double arcsine
    transformation) - if \code{byvar} is not missing.}
  \item{n.w}{Number of observations in subgroups - if \code{byvar} is
    not missing.}
  \item{k.w}{Number of studies combined within subgroups - if
    \code{byvar} is not missing.}
  \item{k.all.w}{Number of all studies in subgroups - if \code{byvar}
    is not missing.}
  \item{Q.w}{Heterogeneity statistics within subgroups - if
    \code{byvar} is not missing.}
  \item{Q.w.fixed}{Overall within subgroups heterogeneity statistic Q
    (based on fixed effect model) - if \code{byvar} is not missing.}
  \item{Q.w.random}{Overall within subgroups heterogeneity statistic Q
    (based on random effects model) - if \code{byvar} is not missing
    (only calculated if argument \code{tau.common} is TRUE).}
  \item{df.Q.w}{Degrees of freedom for test of overall within
    subgroups heterogeneity - if \code{byvar} is not missing.}
  \item{Q.b.fixed}{Overall between subgroups heterogeneity statistic Q
    (based on fixed effect model) - if \code{byvar} is not missing.}
  \item{Q.b.random}{Overall between subgroups heterogeneity statistic
    Q (based on random effects model) - if \code{byvar} is not
    missing.}
  \item{df.Q.b}{Degrees of freedom for test of overall between
    subgroups heterogeneity - if \code{byvar} is not missing.}
  \item{tau.w}{Square-root of between-study variance within subgroups
    - if \code{byvar} is not missing.}
  \item{C.w}{Scaling factor utilised internally to calculate common
    tau-squared across subgroups - if \code{byvar} is not missing.}
  \item{H.w}{Heterogeneity statistic H within subgroups - if
    \code{byvar} is not missing.}
  \item{lower.H.w, upper.H.w}{Lower and upper confidence limti for
    heterogeneity statistic H within subgroups - if \code{byvar} is
    not missing.}
   \item{I2.w}{Heterogeneity statistic I2 within subgroups - if
    \code{byvar} is not missing.}
  \item{lower.I2.w, upper.I2.w}{Lower and upper confidence limti for
    heterogeneity statistic I2 within subgroups - if \code{byvar} is
    not missing.}
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
 
  DerSimonian R & Laird N (1986),
  Meta-analysis in clinical trials. \emph{Controlled Clinical Trials},
  \bold{7}, 177--188.
 
  Higgins JPT, Thompson SG, Spiegelhalter DJ (2009),
  A re-evaluation of random-effects meta-analysis.
  \emph{Journal of the Royal Statistical Society: Series A},
  \bold{172}, 137--159.
  
  Knapp G & Hartung J (2003),
  Improved Tests for a Random Effects Meta-regression with a Single
  Covariate.
  \emph{Statistics in Medicine}, \bold{22}, 2693--710,
  doi: 10.1002/sim.1482 .

  Paule RC & Mandel J (1982),
  Consensus values and weighting factors.
  \emph{Journal of Research of the National Bureau of Standards},
  \bold{87}, 377--385.

  Viechtbauer W (2010),
  Conducting Meta-Analyses in R with the Metafor Package.
  \emph{Journal of Statistical Software}, \bold{36}, 1--48.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{update.meta}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
m1 <- metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))

#
# Print correlations (back transformed from Fisher's z transformation)
#
m1

#
# Print Fisher's z transformed correlations 
#
print(m1, backtransf=FALSE)

#
# Forest plot with back transformed correlations
#
forest(m1)

#
# Forest plot with Fisher's z transformed correlations
#
forest(m1, backtransf=FALSE)


m2 <- update(m1, sm="cor")
m2
# Identical forest plots (as back transformation is the identity transformation)
# forest(m2)
# forest(m2, backtransf=FALSE)
}

\keyword{htest}
