\name{metagen}

\alias{metagen}

\title{Generic inverse variance meta-analysis}

\description{
  Fixed and random effects meta-analysis based on estimates (e.g. log
  hazard ratios) and their standard errors; inverse variance weighting
  is used for pooling.
}

\usage{
metagen(TE, seTE, studlab,
        data=NULL, subset=NULL, exclude=NULL, sm="",
        level=gs("level"), level.comb=gs("level.comb"),
        comb.fixed=gs("comb.fixed"), comb.random=gs("comb.random"),
        hakn=gs("hakn"),
        method.tau=gs("method.tau"), tau.preset=NULL, TE.tau=NULL,
        tau.common=gs("tau.common"),
        prediction=gs("prediction"), level.predict=gs("level.predict"),
        null.effect=0,
        method.bias=gs("method.bias"),
        n.e=NULL, n.c=NULL,
        backtransf=gs("backtransf"),
        pscale=1, irscale = 1, irunit = "person-years",
        title=gs("title"), complab=gs("complab"), outclab="",
        label.e=gs("label.e"), label.c=gs("label.c"),
        label.left=gs("label.left"), label.right=gs("label.right"),
        byvar, bylab, print.byvar=gs("print.byvar"),
        byseparator = gs("byseparator"),
        keepdata=gs("keepdata"),
        warn=gs("warn"))
}

\arguments{
  \item{TE}{Estimate of treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{exclude}{An optional vector specifying studies to exclude from
    meta-analysis, however, to include in printouts and forest plots.}
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
    \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.}
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
  \item{null.effect}{A numeric value specifying the effect under the
    null hypothesis.}
  \item{n.e}{Number of observations in experimental group.}
  \item{n.c}{Number of observations in control group.}
  \item{hakn}{A logical indicating whether method by Hartung and
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
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and plots. If \code{backtransf=TRUE}
    (default), results for \code{sm="OR"} are printed as odds ratios
    rather than log odds ratios and results for \code{sm="ZCOR"} are
    printed as correlations rather than Fisher's z transformed
    correlations, for example.}
  \item{pscale}{A numeric defining a scaling factor for printing of
    single event probabilities, i.e. if argument \code{sm} is equal to
    \code{"PLOGIT"}, \code{"PLN"}, \code{"PRAW"}, \code{"PAS"}, or
    \code{"PFT"}. See also \code{\link{metaprop}}}
  \item{irscale}{A numeric defining a scaling factor for printing of
    rates, i.e. if argument \code{sm} is equal to \code{"IR"},
    \code{"IRLN"}, \code{"IRS"}, or \code{"IRFT"}.}
  \item{irunit}{A character specifying the time unit used to calculate
    rates, e.g. person-years.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{label.e}{Label for experimental group.}
  \item{label.c}{Label for control group.}
  \item{label.left}{Graph label on left side of forest plot.}
  \item{label.right}{Graph label on right side of forest plot.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{TE}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{byseparator}{A character string defining the separator between
    label and levels of grouping variable.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in meta object.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if studies are excluded from meta-analysis due to zero
    standard errors).}
}

\details{
  Generic method for meta-analysis, only treatment estimates and their
  standard error are needed. The method is useful, e.g., for pooling
  of survival data (using log hazard ratio and standard errors as
  input). The inverse variance method is used for pooling.
  
  For several arguments defaults settings are utilised (assignments
  using \code{\link{gs}} function). These defaults can be changed
  using the \code{\link{settings.meta}} function.
  
  Internally, both fixed effect and random effects models are
  calculated regardless of values choosen for arguments
  \code{comb.fixed} and \code{comb.random}. Accordingly, the estimate
  for the random effects model can be extracted from component
  \code{TE.random} of an object of class \code{"meta"} even if
  argument \code{comb.random=FALSE}. However, all functions in R
  package \bold{meta} will adequately consider the values for
  \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.meta}} will not print results for the random
  effects model if \code{comb.random=FALSE}.

  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metagen object by only specifying
  arguments which should be changed.

  For the random effects, the method by Hartung and Knapp (2003) is
  used to adjust test statistics and confidence intervals if argument
  \code{hakn=TRUE}.
  
  The DerSimonian-Laird estimate (1986) is used in the random effects
  model if \code{method.tau="DL"}. The iterative Paule-Mandel method
  (1982) to estimate the between-study variance is used if argument
  \code{method.tau="PM"}. Internally, R function \code{paulemandel} is
  called which is based on R function mpaule.default from R package
  \bold{metRology} from S.L.R. Ellison <s.ellison at lgc.co.uk>.
  
  If R package \bold{metafor} (Viechtbauer 2010) is installed, the
  following methods to estimate the between-study variance
  \eqn{\tau^2} (argument \code{method.tau}) are also available:
  \itemize{
   \item Restricted maximum-likelihood estimator (\code{method.tau="REML"})
   \item Maximum-likelihood estimator (\code{method.tau="ML"})
   \item Hunter-Schmidt estimator (\code{method.tau="HS"})
   \item Sidik-Jonkman estimator (\code{method.tau="SJ"})
   \item Hedges estimator (\code{method.tau="HE"})
   \item Empirical Bayes estimator (\code{method.tau="EB"}).
  }
  For these methods the R function \code{rma.uni} of R package
  \bold{metafor} is called internally. See help page of R function
  \code{rma.uni} for more details on these methods to estimate
  between-study variance.
  
  Argument \code{pscale} can be used to rescale proportions,
  e.g. \code{pscale=1000} means that proportions are expressed as
  events per 1000 observations. This is useful in situations with
  (very) low event probabilities.

  Argument \code{irscale} can be used to rescale rates,
  e.g. \code{irscale=1000} means that rates are expressed as events
  per 1000 time units, e.g. person-years. This is useful in situations
  with (very) low rates. Argument \code{irunit} can be used to specify
  the time unit used in individual studies (default:
  "person-years"). This information is printed in summaries and forest
  plots if argument \code{irscale} is not equal to 1.
}

\value{
  An object of class \code{c("metagen", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{TE, seTE, studlab, exclude, n.e, n.c}{}
  \item{sm, level, level.comb,}{}
  \item{comb.fixed, comb.random,}{}
  \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{}
  \item{tau.common, title, complab, outclab,}{}
  \item{label.e, label.c, label.left, label.right,}{}
  \item{byvar, bylab, print.byvar, byseparator, warn}{As defined above.}
  \item{lower, upper}{Lower and upper confidence interval limits
    for individual studies.}
  \item{zval, pval}{z-value and p-value for test of treatment
    effect for individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
    standard error (fixed effect model).}
  \item{lower.fixed, upper.fixed}{Lower and upper confidence interval limits
  (fixed effect model).}
  \item{zval.fixed, pval.fixed}{z-value and p-value for test of
    overall treatment effect (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall treatment effect and
    standard error (random effects model).}
  \item{lower.random, upper.random}{Lower and upper confidence interval limits
    (random effects model).}
  \item{zval.random, pval.random}{z-value or t-value and corresponding
      p-value for test of overall treatment effect (random effects
      model).}
  \item{prediction, level.predict}{As defined above.}
  \item{seTE.predict}{Standard error utilised for prediction interval.}
  \item{lower.predict, upper.predict}{Lower and upper limits of prediction interval.}
  \item{null.effect}{As defined above.}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic.}
  \item{df.Q}{Degrees of freedom for heterogeneity statistic.}
  \item{tau}{Square-root of between-study variance.}
  \item{se.tau}{Standard error of square-root of between-study variance.}
  \item{C}{Scaling factor utilised internally to calculate common
    tau-squared across subgroups.}
  \item{method}{Pooling method: \code{"Inverse"}.}
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
  \item{n.e.w}{Number of observations in experimental group in
    subgroups - if \code{byvar} is not missing.}
  \item{n.c.w}{Number of observations in control group in subgroups -
    if \code{byvar} is not missing.}
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
  \item{version}{Version of R package \bold{meta} used to create object.}
}

\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.

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
  \emph{Statistics in Medicine}, \bold{22}, 2693--2710,
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

\seealso{\code{\link{update.meta}}, \code{\link{metabin}}, \code{\link{metacont}}, \code{\link{print.meta}}}

\examples{
data(Fleiss93)
meta1 <- metabin(event.e, n.e, event.c, n.c, data=Fleiss93, sm="RR", method="I")
meta1

#
# Identical results by using the following commands:
#
meta1
metagen(meta1$TE, meta1$seTE, sm="RR")

forest(metagen(meta1$TE, meta1$seTE, sm="RR"))


#
# Meta-analysis with prespecified between-study variance
#
summary(metagen(meta1$TE, meta1$seTE, sm="RR", tau.preset=sqrt(0.1)))


#
# Meta-analysis of survival data:
#
logHR <- log(c(0.95, 1.5))
selogHR <- c(0.25, 0.35)

metagen(logHR, selogHR, sm="HR")


#
# Paule-Mandel method to estimate between-study variance
# Data from Paule & Mandel (1982)
#
average <- c(27.044, 26.022, 26.340, 26.787, 26.796)
variance <- c(0.003, 0.076, 0.464, 0.003, 0.014)
#
summary(metagen(average, sqrt(variance), sm="MD", method.tau="PM"))
}

\keyword{htest}
