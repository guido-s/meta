\name{metaprop}

\alias{metaprop}

\title{Meta-analysis of single proportions}

\description{
  Calculation of an overall proportion from studies reporting a single
  proportion. Inverse variance method and generalised linear mixed
  model (GLMM) are available for pooling. For GLMMs, the
  \code{\link[metafor]{rma.glmm}} function from R package
  \bold{metafor} (Viechtbauer 2010) is called internally.
}

\usage{
metaprop(event, n, studlab,
         data=NULL, subset=NULL, method = "Inverse",
         sm=gs("smprop"),
         incr=gs("incr"), allincr=gs("allincr"),
         addincr=gs("addincr"),
         method.ci=gs("method.ci"),
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
         backtransf=gs("backtransf"),
         pscale=1,
         title=gs("title"), complab=gs("complab"), outclab="",
         byvar, bylab, print.byvar=gs("print.byvar"),
         byseparator = gs("byseparator"),
         keepdata=gs("keepdata"),
         warn=gs("warn"),
	 ...)
}

\arguments{
  \item{event}{Number of events.}
  \item{n}{Number of observations.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., event and n.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{method}{A character string indicating which method is to be
    used for pooling of studies. One of \code{"Inverse"} and
    \code{"GLMM"}, can be abbreviated.}
  \item{sm}{A character string indicating which summary measure
    (\code{"PFT"}, \code{"PAS"}, \code{"PRAW"}, \code{"PLN"}, or
    \code{"PLOGIT"}) is to be used for pooling of studies, see Details.}
  \item{incr}{A numeric which is added to event number and sample size
    of studies with zero or all events, i.e., studies with an event
    probability of either 0 or 1.}
  \item{allincr}{A logical indicating if \code{incr} is considered for
    all studies if at least one study has either zero or all
    events. If FALSE (default), \code{incr} is considered only in
    studies with zero or all events.}
  \item{addincr}{A logical indicating if \code{incr} is used for all
    studies irrespective of number of events.}
  \item{method.ci}{A character string indicating which method is used
    to calculate confidence intervals for individual studies, see
    Details.}
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
  \item{method.bias}{A character string indicating which test is to be
    used. Either \code{"rank"}, \code{"linreg"}, or \code{"mm"}, can
    be abbreviated. See function \code{\link{metabias}}.}
  \item{backtransf}{A logical indicating whether results for
    transformed proportions (argument \code{sm!="PRAW"}) should be
    back transformed in printouts and plots. If TRUE (default),
    results will be presented as proportions; otherwise transformed
    proportions will be shown. See Details for presentation of
    confidence intervals.}
  \item{pscale}{A numeric defining a scaling factor for printing of
    single event probabilities.}
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
    to studies with zero or all events should result in a warning.}
  \item{\dots}{Additional arguments passed on to
    \code{\link[metafor]{rma.glmm}} function.}
}

\details{
  Fixed effect and random effects meta-analysis of single proportions
  to calculate an overall proportion. The following transformations of
  proportions are implemented to calculate an overall proportion:
  \itemize{
  \item Logit transformation (\code{sm="PLOGIT"}, default)
  \item Log transformation (\code{sm="PLN"})
  \item Freeman-Tukey Double arcsine transformation (\code{sm="PFT"})  
  \item Arcsine transformation (\code{sm="PAS"})
  \item Raw, i.e. untransformed, proportions (\code{sm="PRAW"})
  }
     
  Note, you should use R function \code{\link{metabin}} to compare
  proportions of pairwise comparisons instead of using \code{metaprop}
  for each treatment arm separately which will break randomisation in
  randomised controlled trials.
  
  Various methods are available to calculate confidence intervals for
  individual study results (see Agresti & Coull 1998; Newcombe 1988):

  \itemize{
  \item Clopper-Pearson interval also called 'exact' binomial interval
    (\code{method.ci="CP"}, default)
  
  \item Wilson Score interval (\code{method.ci="WS"})
  
  \item Wilson Score interval with continuity correction
    (\code{method.ci="WSCC"})
  
  \item Agresti-Coull interval (\code{method.ci="AC"})
  
  \item Simple approximation interval (\code{method.ci="SA"})
  
  \item Simple approximation interval with continuity correction
    (\code{method.ci="SACC"})
  
  \item Normal approximation interval based on summary measure,
    i.e. defined by argument \code{sm} (\code{method.ci="NAsm"})}

  Note, with exception of the normal approximation based on the
  summary measure, i.e. \code{method.ci="NAsm"}, the same confidence
  interval is calculated for any summary measure (argument \code{sm})
  as only number of events and observations are used in the
  calculation disregarding the chosen summary measure. Results will be
  presented for transformed proportions if argument
  \code{backtransf=FALSE} in the \code{\link{print.meta}},
  \code{\link{print.summary.meta}}, or \code{\link{forest.meta}}
  function. In this case, argument \code{method.ci="NAsm"} is used,
  i.e. confidence intervals based on the normal approximation based on
  the summary measure.
  
  Argument \code{pscale} can be used to rescale proportions,
  e.g. \code{pscale=1000} means that proportions are expressed as
  events per 1000 observations. This is useful in situations with
  (very) low event probabilities.
  
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
  
  A distinctive and frequently overlooked advantage of binary data is
  that individual patient data (IPD) can be extracted. Accordingly, a
  random intercept logistic regression model can be utilised for the
  meta-analysis of proportions (Stijnen et al., 2010). This method is
  available (argument \code{method = "GLMM"}) by calling the
  \code{\link[metafor]{rma.glmm}} function from R package
  \bold{metafor} internally.

  If the summary measure is equal to "PRAW", "PLN", or "PLOGIT", a
  continuity correction is applied if any study has either zero or all
  events, i.e., an event probability of either 0 or 1. By default, 0.5
  is used as continuity correction (argument \code{incr}). This
  continuity correction is used both to calculate individual study
  results with confidence limits and to conduct meta-analysis based on
  the inverse variance method. For GLMMs no continuity correction is
  used.
  
  Argument \code{byvar} can be used to conduct subgroup analysis for
  all methods but GLMMs. Instead use the \code{\link{metareg}}
  function for GLMMs which can also be used for continuous covariates.

  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metaprop object by only specifying
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
  An object of class \code{c("metaprop", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{event, n, studlab,}{}
  \item{sm, incr, allincr, addincr, method.ci,}{}
  \item{level, level.comb,}{As defined above.}
  \item{comb.fixed, comb.random,}{}
  \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{}
  \item{tau.common, title, complab, outclab,}{}
  \item{byvar, bylab, print.byvar, byseparator, warn}{}
  \item{TE, seTE}{Estimated (un)transformed proportion and its standard
    error for individual studies.}
  \item{lower, upper}{Lower and upper confidence interval limits
    for individual studies.}
  \item{zval, pval}{z-value and p-value for test of treatment
    effect for individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall (un)transformed
    proportion and standard error (fixed effect model).}
  \item{lower.fixed, upper.fixed}{Lower and upper confidence interval limits
  (fixed effect model).}
  \item{zval.fixed, pval.fixed}{z-value and p-value for test of
    overall effect (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall (un)transformed
    proportion and standard error (random effects model).}
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
  Agresti A & Coull BA (1998),
  Approximate is better than “exact” for interval estimation of
  binomial proportions.
  \emph{The American Statistician},
  \bold{52}, 119--126.

  DerSimonian R & Laird N (1986),
  Meta-analysis in clinical trials.
  \emph{Controlled Clinical Trials},
  \bold{7}, 177--188.

  Edward JM et al. (2006),
  Adherence to antiretroviral therapy in sub-saharan Africa and North
  America - a meta-analysis.
  \emph{Journal of the American Medical Association},
  \bold{296}, 679--690.

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

  Miller JJ (1978),
  The inverse of the Freeman-Tukey double arcsine transformation.
  \emph{The American Statistician}, \bold{32}, 138.
  
  Newcombe RG (1998),
  Two-sided confidence intervals for the single proportion: Comparison
  of seven methods.
  \emph{Statistics in Medicine}, \bold{17}, 857--872.
  
  Paule RC & Mandel J (1982),
  Consensus values and weighting factors.
  \emph{Journal of Research of the National Bureau of Standards},
  \bold{87}, 377--385.

  Pettigrew HM, Gart JJ, Thomas DG (1986),
  The bias and higher cumulants of the logarithm of a binomial
  variate. \emph{Biometrika}, \bold{73}, 425--435.
  
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
# Apply various meta-analysis methods to estimate proportions
#
m1 <- metaprop(4:1, c(10, 20, 30, 40))
m2 <- update(m1, sm="PAS")
m3 <- update(m1, sm="PRAW")
m4 <- update(m1, sm="PLN")
m5 <- update(m1, sm="PFT")
#
m1
m2
m3
m4
m5
#
forest(m1)
# forest(m2)
# forest(m3)
# forest(m3, pscale=100)
# forest(m4)
# forest(m5)


#
# Do not back transform results, e.g. print logit transformed
# proportions if sm="PLOGIT" and store old settings
#
oldset <- settings.meta(backtransf=FALSE)
#
m6  <- metaprop(4:1, c(10, 20, 30, 40))
m7  <- update(m6, sm="PAS")
m8  <- update(m6, sm="PRAW")
m9  <- update(m6, sm="PLN")
m10 <- update(m6, sm="PFT")
#
forest(m6)
# forest(m7)
# forest(m8)
# forest(m8, pscale=100)
# forest(m9)
# forest(m10)
#
# Use old settings
#
settings.meta(oldset)


#
# Examples with zero events
#
m1 <- metaprop(c(0, 0, 10, 10), rep(100, 4))
m2 <- metaprop(c(0, 0, 10, 10), rep(100, 4), incr=0.1)
#
summary(m1)
summary(m2)
#
# forest(m1)
# forest(m2)


#
# Example from Miller (1978):
#
death <- c(3, 6, 10, 1)
animals <- c(11, 17, 21, 6)
#
m3 <- metaprop(death, animals, sm="PFT")
forest(m3)


#
# Data examples from Newcombe (1998)
# - apply various methods to estimate confidence intervals for
#   individual studies
#
event <- c(81, 15, 0, 1)
n <- c(263, 148, 20, 29)
#
m1 <- metaprop(event, n, sm="PLOGIT", method.ci="SA")
m2 <- update(m1, method.ci="SACC")
m3 <- update(m1, method.ci="WS")
m4 <- update(m1, method.ci="WSCC")
m5 <- update(m1, method.ci="CP")
#
lower <- round(rbind(NA, m1$lower, m2$lower, NA, m3$lower, m4$lower, NA, m5$lower), 4)
upper <- round(rbind(NA, m1$upper, m2$upper, NA, m3$upper, m4$upper, NA, m5$upper), 4)
#
tab1 <- data.frame(
  scen1=meta:::p.ci(lower[,1], upper[,1]),
  scen2=meta:::p.ci(lower[,2], upper[,2]),
  scen3=meta:::p.ci(lower[,3], upper[,3]),
  scen4=meta:::p.ci(lower[,4], upper[,4]),
  stringsAsFactors=FALSE
  )
names(tab1) <- c("r=81, n=263", "r=15, n=148", "r=0, n=20", "r=1, n=29")
row.names(tab1) <- c("Simple", "- SA", "- SACC",
                     "Score", "- WS", "- WSCC",
                     "Binomial", "- CP")
tab1[is.na(tab1)] <- ""
#
# Newcombe (1998), Table I, methods 1-5:
#
tab1
#
# Same confidence interval, i.e. unaffected by choice of summary measure
#
print(metaprop(event, n, sm="PLOGIT", method.ci="WS"), ma=FALSE)
print(metaprop(event, n, sm="PLN", method.ci="WS"), ma=FALSE)
print(metaprop(event, n, sm="PFT", method.ci="WS"), ma=FALSE)
print(metaprop(event, n, sm="PAS", method.ci="WS"), ma=FALSE)
print(metaprop(event, n, sm="PRAW", method.ci="WS"), ma=FALSE)
#
# Different confidence intervals as argument sm="NAsm"
#
print(metaprop(event, n, sm="PLOGIT", method.ci="NAsm"), ma=FALSE)
print(metaprop(event, n, sm="PLN", method.ci="NAsm"), ma=FALSE)
print(metaprop(event, n, sm="PFT", method.ci="NAsm"), ma=FALSE)
print(metaprop(event, n, sm="PAS", method.ci="NAsm"), ma=FALSE)
print(metaprop(event, n, sm="PRAW", method.ci="NAsm"), ma=FALSE)
#
# Different confidence intervals as argument backtransf=FALSE.
# Accordingly, method.ci="NAsm" used internally.
#
print(metaprop(event, n, sm="PLOGIT", method.ci="WS"), ma=FALSE, backtransf=FALSE)
print(metaprop(event, n, sm="PLN", method.ci="WS"), ma=FALSE, backtransf=FALSE)
print(metaprop(event, n, sm="PFT", method.ci="WS"), ma=FALSE, backtransf=FALSE)
print(metaprop(event, n, sm="PAS", method.ci="WS"), ma=FALSE, backtransf=FALSE)
print(metaprop(event, n, sm="PRAW", method.ci="WS"), ma=FALSE, backtransf=FALSE)
#
# Same results (printed on original and log scale, respectively)
#
print(metaprop(event, n, sm="PLN", method.ci="NAsm"), ma=FALSE)
print(metaprop(event, n, sm="PLN"), ma=FALSE, backtransf=FALSE)
# Results for first study (on log scale)
round(log(c(0.3079848, 0.2569522, 0.3691529)), 4)


#
# Meta-analysis using generalised linear mixed models
# (only if R packages 'metafor' and 'lme4' are available)
#
if (suppressMessages(require(metafor, quietly = TRUE, warn = FALSE)) &
    require(lme4, quietly = TRUE))
 metaprop(event, n, method = "GLMM")


#
# Print results as events per 1000 observations
#
print(metaprop(6:8, c(100, 1200, 1000)), pscale = 1000, digits = 1)
}

\keyword{htest}
