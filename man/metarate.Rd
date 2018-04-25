\name{metarate}

\alias{metarate}

\title{Meta-analysis of single incidence rates}

\description{
  Calculation of an overall incidence rate from studies reporting a
  single incidence rate. Inverse variance method and generalised
  linear mixed model (GLMM) are available for pooling. For GLMMs, the
  \code{\link[metafor]{rma.glmm}} function from R package
  \bold{metafor} (Viechtbauer 2010) is called internally.
}

\usage{
metarate(event, time, studlab,
         data=NULL, subset=NULL, exclude=NULL,
         method = "Inverse",
         sm=gs("smrate"),
         incr=gs("incr"), allincr=gs("allincr"),
         addincr=gs("addincr"),
         level=gs("level"), level.comb=gs("level.comb"),
         comb.fixed=gs("comb.fixed"), comb.random=gs("comb.random"),
         hakn=gs("hakn"),
         method.tau=
         ifelse(!is.na(charmatch(tolower(method), "glmm", nomatch = NA)),
                "ML", gs("method.tau")),
         tau.preset=NULL, TE.tau=NULL,
         tau.common=gs("tau.common"),
         prediction=gs("prediction"), level.predict=gs("level.predict"),
	 null.effect=NA,
         method.bias=gs("method.bias"),
         backtransf=gs("backtransf"),
         irscale=1, irunit="person-years",
         title=gs("title"), complab=gs("complab"), outclab="",
         byvar, bylab, print.byvar=gs("print.byvar"),
         byseparator = gs("byseparator"),
         keepdata=gs("keepdata"),
         warn=gs("warn"),
	 ...)
}

\arguments{
  \item{event}{Number of events.}
  \item{time}{Person time at risk.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., event and time.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{exclude}{An optional vector specifying studies to exclude from
    meta-analysis, however, to include in printouts and forest plots.}
  \item{method}{A character string indicating which method is to be
    used for pooling of studies. One of \code{"Inverse"} and
    \code{"GLMM"}, can be abbreviated.}
  \item{sm}{A character string indicating which summary measure
    (\code{"IR"}, \code{"IRLN"}, \code{"IRS"}, or \code{"IRFT"}) is to
    be used for pooling of studies, see Details.}
  \item{incr}{A numeric which is added to the event number of studies
    with zero events, i.e., studies with an incidence rate of 0.}
  \item{allincr}{A logical indicating if \code{incr} is considered for
    all studies if at least one study has zero events. If FALSE
    (default), \code{incr} is considered only in studies with zero
    events.}
  \item{addincr}{A logical indicating if \code{incr} is used for all
    studies irrespective of number of events.}
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
    Knapp should be used to adjust test statistics and confidence
    intervals.}
  \item{method.tau}{A character string indicating which method is used
    to estimate the between-study variance \eqn{\tau^2}, see Details.}
  \item{tau.preset}{Prespecified value for the square-root of the
    between-study variance \eqn{\tau^2}.}
  \item{TE.tau}{Overall treatment effect used to estimate the
    between-study variance tau-squared.}
  \item{tau.common}{A logical indicating whether tau-squared should be
    the same across subgroups.}
  \item{null.effect}{A numeric value specifying the effect under the
    null hypothesis.}
  \item{method.bias}{A character string indicating which test is to be
    used. Either \code{"rank"}, \code{"linreg"}, or \code{"mm"}, can
    be abbreviated. See function \code{\link{metabias}}.}
  \item{backtransf}{A logical indicating whether results for
    transformed rates (argument \code{sm!="IR"}) should be back
    transformed in printouts and plots. If TRUE (default), results
    will be presented as incidence rates; otherwise transformed
    rates will be shown.}
  \item{irscale}{A numeric defining a scaling factor for printing of
    rates.}
  \item{irunit}{A character string specifying the time unit used to
    calculate rates, e.g. person-years.}    
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{event}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{byseparator}{A character string defining the separator between
    label and levels of grouping variable.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in meta object.}
  \item{warn}{A logical indicating whether the addition of \code{incr}
    to studies with zero events should result in a warning.}
  \item{\dots}{Additional arguments passed on to
    \code{\link[metafor]{rma.glmm}} function.}
}

\details{
  Fixed effect and random effects meta-analysis of single incidence
  rates to calculate an overall rate. The following transformations of
  incidence rates are implemented to calculate an overall rate:
  \itemize{
  \item Log transformation (\code{sm="IRLN"}, default)
  \item Square root transformation (\code{sm="IRS"})
  \item Freeman-Tukey Double arcsine transformation (\code{sm="IRFT"})  
  \item No transformation (\code{sm="IR"})
  }
  
  Note, you should use R function \code{\link{metainc}} to compare
  incidence rates of pairwise comparisons instead of using
  \code{metarate} for each treatment arm separately which will break
  randomisation in randomised controlled trials.
  
  Argument \code{irscale} can be used to rescale rates,
  e.g. \code{irscale=1000} means that rates are expressed as events
  per 1000 time units, e.g. person-years. This is useful in situations
  with (very) low rates. Argument \code{irunit} can be used to specify
  the time unit used in individual studies (default:
  "person-years"). This information is printed in summaries and forest
  plots if argument \code{irscale} is not equal to 1.
  
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
  
  A random intercept Poisson regression model can be utilised for the
  meta-analysis of incidence rates (Stijnen et al., 2010). This method
  is available (argument \code{method = "GLMM"}) by calling the
  \code{\link[metafor]{rma.glmm}} function from R package
  \bold{metafor} internally.

  If the summary measure is equal to "IR" or "IRLN", a continuity
  correction is applied if any study has zero events, i.e., an
  incidence rate of 0. By default, 0.5 is used as continuity
  correction (argument \code{incr}). This continuity correction is
  used both to calculate individual study results with confidence
  limits and to conduct meta-analysis based on the inverse variance
  method. For Freeman-Tukey and square root transformation and GLMMs
  no continuity correction is used.
  
  Argument \code{byvar} can be used to conduct subgroup analysis for
  all methods but GLMMs. Instead use the \code{\link{metareg}}
  function for GLMMs which can also be used for continuous covariates.

  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metarate object by only specifying
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
}

\value{
  An object of class \code{c("metarate", "meta")} with corresponding
  \code{print}, \code{summary}, and \code{forest} functions. The
  object is a list containing the following components:
  \item{event, n, studlab, exclude,}{}
  \item{sm, incr, allincr, addincr, method.ci,}{}
  \item{level, level.comb,}{As defined above.}
  \item{comb.fixed, comb.random,}{}
  \item{hakn, method.tau, tau.preset, TE.tau, null.effect,}{}
  \item{method.bias, tau.common, title, complab, outclab,}{}
  \item{byvar, bylab, print.byvar, byseparator, warn}{}
  \item{TE, seTE}{Estimated (un)transformed incidence rate and its
    standard error for individual studies.}
  \item{lower, upper}{Lower and upper confidence interval limits
    for individual studies.}
  \item{zval, pval}{z-value and p-value for test of treatment
    effect for individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall (un)transformed
    incidence rate and standard error (fixed effect model).}
  \item{lower.fixed, upper.fixed}{Lower and upper confidence interval limits
  (fixed effect model).}
  \item{zval.fixed, pval.fixed}{z-value and p-value for test of
    overall effect (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall (un)transformed
    incidence rate and standard error (random effects model).}
  \item{lower.random, upper.random}{Lower and upper confidence interval limits
  (random effects model).}
  \item{zval.random, pval.random}{z-value or t-value and corresponding
    p-value for test of overall effect (random effects
    model).}
  \item{prediction, level.predict}{As defined above.}
  \item{seTE.predict}{Standard error utilised for prediction interval.}
  \item{lower.predict, upper.predict}{Lower and upper limits of prediction interval.}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{df.Q}{Degrees of freedom for heterogeneity statistic.}
  \item{pval.Q}{P-value of heterogeneity test.}
  \item{Q.LRT}{Heterogeneity statistic for likelihood-ratio test (only
    if \code{method = "GLMM"}).}
  \item{df.Q.LRT}{Degrees of freedom for likelihood-ratio test}
  \item{pval.Q.LRT}{P-value of likelihood-ratio test.}
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
  \item{event.w}{Number of events in subgroups - if \code{byvar} is
    not missing.}
  \item{n.w}{Number of observations in subgroups - if \code{byvar} is
    not missing.}
  \item{k.w}{Number of studies combined within subgroups - if
    \code{byvar} is not missing.}
  \item{k.all.w}{Number of all studies in subgroups - if \code{byvar}
    is not missing.}
  \item{Q.w.fixed}{Overall within subgroups heterogeneity statistic Q
    (based on fixed effect model) - if \code{byvar} is not missing.}    
  \item{Q.w.random}{Overall within subgroups heterogeneity statistic Q
    (based on random effects model) - if \code{byvar} is not missing
    (only calculated if argument \code{tau.common} is TRUE).}
  \item{df.Q.w}{Degrees of freedom for test of overall within
    subgroups heterogeneity - if \code{byvar} is not missing.}
  \item{pval.Q.w.fixed}{P-value of within subgroups heterogeneity
    statistic Q (based on fixed effect model) - if \code{byvar} is not
    missing.}
  \item{pval.Q.w.random}{P-value of within subgroups heterogeneity
    statistic Q (based on random effects model) - if \code{byvar} is
    not missing.}
  \item{Q.b.fixed}{Overall between subgroups heterogeneity statistic Q
    (based on fixed effect model) - if \code{byvar} is not missing.}
  \item{Q.b.random}{Overall between subgroups heterogeneity statistic
    Q (based on random effects model) - if \code{byvar} is not
    missing.}
  \item{df.Q.b}{Degrees of freedom for test of overall between
    subgroups heterogeneity - if \code{byvar} is not missing.}
  \item{pval.Q.b.fixed}{P-value of between subgroups heterogeneity
    statistic Q (based on fixed effect model) - if \code{byvar} is not
    missing.}
  \item{pval.Q.b.random}{P-value of between subgroups heterogeneity
    statistic Q (based on random effects model) - if \code{byvar} is
    not missing.}
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
  \item{incr.event}{Increment added to number of events.}
  \item{keepdata}{As defined above.}
  \item{data}{Original data (set) used in function call (if
    \code{keepdata=TRUE}).}
  \item{subset}{Information on subset of original data used in
    meta-analysis (if \code{keepdata=TRUE}).}
  \item{.glmm.fixed}{GLMM object generated by call of
    \code{\link[metafor]{rma.glmm}} function (fixed effect model).}
  \item{.glmm.random}{GLMM object generated by call of
    \code{\link[metafor]{rma.glmm}} function (random effects model).}
  \item{call}{Function call.}
  \item{version}{Version of R package \bold{meta} used to create object.}
  \item{version.metafor}{Version of R package \bold{metafor} used for GLMMs.}
}

\references{
  DerSimonian R & Laird N (1986),
  Meta-analysis in clinical trials.
  \emph{Controlled Clinical Trials},
  \bold{7}, 177--188.

  Freeman MF & Tukey JW (1950),
  Transformations related to the angular and the square
  root. \emph{Annals of Mathematical Statistics},
  \bold{21}, 607--611.
  
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
  
  Stijnen T, Hamza TH, Ozdemir P (2010),
  Random effects meta-analysis of event outcome in the framework of
  the generalized linear mixed model with applications in sparse
  data.
  \emph{Statistics in Medicine}, \bold{29}, 3046--67.
  
  Viechtbauer W (2010),
  Conducting Meta-Analyses in R with the Metafor Package.
  \emph{Journal of Statistical Software}, \bold{36}, 1--48.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{update.meta}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
#
# Apply various meta-analysis methods to estimate incidence rates
#
m1 <- metarate(4:1, c(10, 20, 30, 40))
m2 <- update(m1, sm="IR")
m3 <- update(m1, sm="IRS")
m4 <- update(m1, sm="IRFT")
#
m1
m2
m3
m4
#
forest(m1)
forest(m1, irscale=100)
forest(m1, irscale=100, irunit="person-days")
forest(m1, backtransf = FALSE)
# forest(m2)
# forest(m3)
# forest(m4)
#
m5 <- metarate(40:37, c(100, 200, 300, 400), sm="IRFT")
m5
}

\keyword{htest}
