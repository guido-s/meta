\name{metainc}

\alias{metainc}

\title{Meta-analysis of incidence rates}

\description{
  Calculation of fixed effect and random effects estimates (incidence
  rate ratio or incidence rate difference) for meta-analyses with
  event counts. Mantel-Haenszel, Cochran, inverse variance method, and
  generalised linear mixed model (GLMM) are available for pooling. For
  GLMMs, the \code{\link[metafor]{rma.glmm}} function from R package
  \bold{metafor} (Viechtbauer 2010) is called internally.
}

\usage{
metainc(event.e, time.e, event.c, time.c, studlab,
        data=NULL, subset=NULL, method="MH",
        sm=gs("sminc"),
        incr=gs("incr"), allincr=gs("allincr"),
        addincr=gs("addincr"),
        model.glmm = "UM.FS",
        level=gs("level"), level.comb=gs("level.comb"),
        comb.fixed=gs("comb.fixed"), comb.random=gs("comb.random"),
        hakn=gs("hakn"),
        method.tau=
        ifelse(!is.na(charmatch(tolower(method), "glmm", nomatch = NA)),
               "ML", gs("method.tau")),
        tau.preset=NULL, TE.tau=NULL,
        tau.common=gs("tau.common"),
        prediction=gs("prediction"), level.predict=gs("level.predict"),
        method.bias=gs("method.bias"),
        n.e=NULL, n.c=NULL,
        backtransf=gs("backtransf"),
        title=gs("title"), complab=gs("complab"), outclab="",
        label.e=gs("label.e"), label.c=gs("label.c"),
        label.left=gs("label.left"), label.right=gs("label.right"),
        byvar, bylab, print.byvar=gs("print.byvar"),
        byseparator = gs("byseparator"),
        keepdata=gs("keepdata"),
        warn=gs("warn"),
	...)
}

\arguments{
  \item{event.e}{Number of events in experimental group.}
  \item{time.e}{Person time at risk in experimental group.}
  \item{event.c}{Number of events in control group.}
  \item{time.c}{Person time at risk in control group.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., event.e, time.e, event.c, and time.c.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{method}{A character string indicating which method is to be
    used for pooling of studies. One of \code{"MH"}, \code{"Inverse"},
    \code{"Cochran"}, or \code{"GLMM"} can be abbreviated.}
  \item{sm}{A character string indicating which summary measure
    (\code{"IRR"} or \code{"IRD"}) is to be used for pooling of
    studies, see Details.}
  \item{incr}{A numerical value which is added to each cell frequency
    for studies with a zero cell count, see Details.}
  \item{allincr}{A logical indicating if \code{incr} is added to each
    cell frequency of all studies if at least one study has a zero cell
    count. If FALSE (default), \code{incr} is added only to each cell frequency of
    studies with a zero cell count.}
  \item{addincr}{A logical indicating if \code{incr} is added to each cell
    frequency of all studies irrespective of zero cell counts.}
  \item{model.glmm}{A character string indicating which GLMM should be
    used. One of \code{"UM.FS"}, \code{"UM.RS"}, and \code{"CM.EL"},
    see Details.}
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
  between-study variance \eqn{\tau^2}.}
  \item{tau.common}{A logical indicating whether tau-squared should be
    the same across subgroups.}
  \item{method.bias}{A character string indicating which test for
    funnel plot asymmetry is to be used. Either \code{"linreg"} or
    \code{"rank"}, can be abbreviated. See function \code{\link{metabias}}}
  \item{n.e}{Number of observations in experimental group (optional).}
  \item{n.c}{Number of observations in control group (optional).}
  \item{backtransf}{A logical indicating whether results for incidence
    rate ratio (\code{sm="IRR"}) should be back transformed in
    printouts and plots. If TRUE (default), results will be presented
    as incidence rate ratios; otherwise log incidence rate ratios will
    be shown.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{label.e}{Label for experimental group.}
  \item{label.c}{Label for control group.}
  \item{label.left}{Graph label on left side of forest plot.}
  \item{label.right}{Graph label on right side of forest plot.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{event.e}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{byseparator}{A character string defining the separator between
    label and levels of grouping variable.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in meta object.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if \code{incr} is added to studies with zero cell
    frequencies).}
  \item{\dots}{Additional arguments passed on to
    \code{\link[metafor]{rma.glmm}} function.}
}

\details{
  Treatment estimates and standard errors are calculated for each
  study. The following measures of treatment effect are available:
  \itemize{
   \item Incidence Rate Ratio (\code{sm="IRR"})
   \item Incidence Rate Difference (\code{sm="IRD"})
  }
    
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

  By default, both fixed effect and random effects models are
  considered (see arguments \code{comb.fixed} and
  \code{comb.random}). If \code{method} is \code{"MH"} (default), the
  Mantel-Haenszel method is used to calculate the fixed effect
  estimate (Greenland & Robbins, 1985); if \code{method} is
  \code{"Inverse"}, inverse variance weighting is used for pooling; if
  \code{method} is \code{"Cochran"}, the Cochran method is used for
  pooling (Bayne-Jones, 1964, Chapter 8).
  
  A distinctive and frequently overlooked advantage of incidence rates
  is that individual patient data (IPD) can be extracted from count
  data. Accordingly, statistical methods for IPD, i.e., generalised
  linear mixed models, can be utilised in a meta-analysis of incidence
  rate ratios (Stijnen et al., 2010). These methods are available
  (argument \code{method = "GLMM"}) by calling the
  \code{\link[metafor]{rma.glmm}} function from R package
  \bold{metafor} internally. Three different GLMMs are available for
  meta-analysis of incidence rate ratios using argument
  \code{model.glmm} (which corresponds to argument \code{model} in the
  \code{\link[metafor]{rma.glmm}} function):
  \itemize{
    \item Poisson regression model with fixed study effects (default)
    \item[] (\code{model.glmm = "UM.FS"}, i.e., \bold{U}nconditional
      \bold{M}odel - \bold{F}ixed \bold{S}tudy effects)
    \item Mixed-effects Poisson regression model with random study
      effects
    \item[] (\code{model.glmm = "UM.RS"}, i.e., \bold{U}nconditional
      \bold{M}odel - \bold{R}andom \bold{S}tudy effects)
    \item Generalised linear mixed model (conditional Poisson-Normal)
    \item[] (\code{model.glmm = "CM.EL"}, i.e., \bold{C}onditional
      \bold{M}odel - \bold{E}xact \bold{L}ikelihood)
  }
  Details on these three GLMMs as well as additional arguments which
  can be provided using argument '\code{\dots}' in \code{metainc} are
  described in \code{\link[metafor]{rma.glmm}} where you can also find
  information on the iterative algorithms used for estimation.  Note,
  regardless of which value is used for argument \code{model.glmm},
  results for two different GLMMs are calculated: fixed effect model
  (with fixed treatment effect) and random effects model (with random
  treatment effects).
  
  For studies with a zero cell count, by default, 0.5 is added to all
  cell frequencies of these studies (argument \code{incr}). This
  continuity correction is used both to calculate individual study
  results with confidence limits and to conduct meta-analysis based on
  the inverse variance method. For Mantel-Haenszel method, Cochran
  method, and GLMMs, nothing is added to zero cell counts.
  Accordingly, estimates for these methods are not defined if the
  number of events is zero in all studies either in the experimental
  or control group.
  
  Argument \code{byvar} can be used to conduct subgroup analysis for
  all methods but GLMMs. Instead use the \code{\link{metareg}}
  function for GLMMs which can also be used for continuous covariates.
  
  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metainc object by only specifying
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
  An object of class \code{c("metainc", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{event.e, time.e, event.c, time.c, studlab,}{}
  \item{sm, method, incr, allincr, addincr, model.glmm, warn,}{}
  \item{level, level.comb, comb.fixed, comb.random,}{}
  \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{}
  \item{tau.common, title, complab, outclab,}{}
  \item{label.e, label.c, label.left, label.right,}{}
  \item{byvar, bylab, print.byvar, byseparator}{As defined above.}
  \item{TE, seTE}{Estimated treatment effect and standard error of individual studies.}
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
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{df.Q}{Degrees of freedom for heterogeneity statistic.}
  \item{tau}{Square-root of between-study variance.}
  \item{se.tau}{Standard error of square-root of between-study variance.}
  \item{C}{Scaling factor utilised internally to calculate common
    tau-squared across subgroups.}
  \item{sparse}{Logical flag indicating if any study included in
    meta-analysis has any zero cell frequencies.}
  \item{incr.event}{Increment added to number of events.}
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
  \item{event.e.w}{Number of events in experimental group in subgroups
    - if \code{byvar} is not missing.}
  \item{time.e.w}{Total person time in subgroups (experimental group)
    - if \code{byvar} is not missing.}
  \item{n.e.w}{Number of observations in experimental group in
    subgroups - if \code{byvar} is not missing.}
  \item{event.c.w}{Number of events in control group in subgroups - if
    \code{byvar} is not missing.}
  \item{time.c.w}{Total person time in subgroups (control group) - if
    \code{byvar} is not missing.}
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
  \item{.glmm.fixed}{GLMM object generated by call of
    \code{\link[metafor]{rma.glmm}} function (fixed effect model).}
  \item{.glmm.random}{GLMM object generated by call of
    \code{\link[metafor]{rma.glmm}} function (random effects model).}
  \item{call}{Function call.}
  \item{version}{Version of R package \bold{meta} used to create object.}
  \item{version.metafor}{Version of R package \bold{metafor} used for GLMMs.}
}

\references{
  Bayne-Jones S et al. (1964),
  Smoking and Health: Report of the Advisory Committee to the Surgeon
  General of the United States. U-23 Department of Health, Education,
  and Welfare. Public Health Service Publication No. 1103.
  \url{http://profiles.nlm.nih.gov/ps/retrieve/ResourceMetadata/NNBBMQ}

  DerSimonian R & Laird N (1986),
  Meta-analysis in clinical trials. \emph{Controlled Clinical Trials},
  \bold{7}, 177--188.

  Greenland S & Robins JM (1985),
  Estimation of a common effect parameter from sparse follow-up data.
  \emph{Biometrics}, \bold{41}, 55--68.

  Hartung J & Knapp G (2001),
  A Refined Method for the Meta-analysis of Controlled Clinical Trials
  with Binary Outcome.
  \emph{Statistics in Medicine}, \bold{20}, 3875--89.
  
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

\seealso{\code{\link{metabin}}, \code{\link{update.meta}}, \code{\link{print.meta}}}

\examples{
data(smoking)

m1 <- metainc(d.smokers, py.smokers,
              d.nonsmokers, py.nonsmokers,
              data=smoking, studlab=study)
print(m1, digits=2)

m2 <- metainc(d.smokers, py.smokers,
              d.nonsmokers, py.nonsmokers,
              data=smoking, studlab=study,
              method="Cochran")
print(m2, digits=2)

data(lungcancer)

m3 <- metainc(d.smokers, py.smokers,
              d.nonsmokers, py.nonsmokers,
              data=lungcancer, studlab=study)
print(m3, digits=2)


# Redo Cochran meta-analysis with inflated standard errors
#
# All cause mortality
#
TEa <- log( (smoking$d.smokers/smoking$py.smokers) /
            (smoking$d.nonsmokers/smoking$py.nonsmokers)
          )
seTEa <- sqrt(1/smoking$d.smokers +
              1/smoking$d.nonsmokers + 2.5/smoking$d.nonsmokers)
#
metagen(TEa, seTEa, sm="IRR", studlab=smoking$study)

# Lung cancer mortality
#
TEl <- log( (lungcancer$d.smokers/lungcancer$py.smokers) /
            (lungcancer$d.nonsmokers/lungcancer$py.nonsmokers)
          )
seTEl <- sqrt(1/lungcancer$d.smokers +
              1/lungcancer$d.nonsmokers + 2.25/lungcancer$d.nonsmokers)
#
metagen(TEl, seTEl, sm="IRR", studlab=lungcancer$study)


#
# Meta-analysis using generalised linear mixed models
# (only if R packages 'metafor' and 'lme4' are available)
#
if (suppressMessages(require(metafor, quietly = TRUE, warn = FALSE)) &
    require(lme4, quietly = TRUE))
 #
 # Poisson regression model (fixed study effects)
 #
 m4 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
               data = smoking, studlab = study, method = "GLMM")
 m4
 #
 # Mixed-effects Poisson regression model (random study effects)
 #
 update(m4, model.glmm = "UM.RS", nAGQ = 1)
 #
 # Generalised linear mixed model (conditional Poisson-Normal)
 #
 update(m4, model.glmm = "CM.EL")
}
\keyword{htest}
