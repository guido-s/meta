\name{metabin}

\alias{metabin}

\title{Meta-analysis of binary outcome data}

\description{
  Calculation of fixed and random effects estimates (risk ratio, odds
  ratio, risk difference or arcsine difference) for meta-analyses with
  binary outcome data. Mantel-Haenszel, inverse variance and Peto
  method are available for pooling.
}

\usage{
metabin(event.e, n.e, event.c, n.c, studlab,
        data=NULL, subset=NULL,
        method=ifelse(tau.common, "Inverse", .settings$method),
        sm=
        ifelse(!is.na(charmatch(method, c("Peto", "peto"),
                                nomatch = NA)),
               "OR", .settings$smbin),
        incr=.settings$incr, allincr=.settings$allincr,
        addincr=.settings$addincr, allstudies=.settings$allstudies,
        MH.exact=.settings$MH.exact, RR.cochrane=.settings$RR.cochrane,
        level=.settings$level, level.comb=.settings$level.comb,
        comb.fixed=.settings$comb.fixed, comb.random=.settings$comb.random,
        hakn=.settings$hakn,
        method.tau=.settings$method.tau, tau.preset=NULL, TE.tau=NULL,
        tau.common=.settings$tau.common,
        prediction=.settings$prediction, level.predict=.settings$level.predict,
        method.bias=ifelse(sm=="OR", "score", .settings$method.bias),
        backtransf=.settings$backtransf,
        title=.settings$title, complab=.settings$complab, outclab="",
        label.e=.settings$label.e, label.c=.settings$label.c,
        label.left=.settings$label.left, label.right=.settings$label.right,
        byvar, bylab, print.byvar=.settings$print.byvar,
        print.CMH=.settings$print.CMH,
        keepdata=.settings$keepdata,
        warn=.settings$warn)
}

\arguments{
  \item{event.e}{Number of events in experimental group.}
  \item{n.e}{Number of observations in experimental group.}
  \item{event.c}{Number of events in control group.}
  \item{n.c}{Number of observations in control group.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., event.e, n.e, event.c, and n.c.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{method}{A character string indicating which method is to be used
    for pooling of studies. One of \code{"Inverse"}, \code{"MH"}, or
    \code{"Peto"}, can be abbreviated.}
  \item{sm}{A character string indicating which summary measure
    (\code{"RR"}, \code{"OR"}, \code{"RD"}, or \code{"ASD"}) is to be used
    for pooling of studies, see Details.}
  \item{incr}{Could be either a numerical value which is added to each
    cell frequency for studies with a zero cell count or the character
    string \code{"TACC"} which stands for treatment arm continuity
    correction, see Details.}
  \item{allincr}{A logical indicating if \code{incr} is added to each
    cell frequency of all studies if at least one study has a zero cell
    count. If FALSE (default), \code{incr} is added only to each cell frequency of
    studies with a zero cell count.}
  \item{addincr}{A logical indicating if \code{incr} is added to each cell
    frequency of all studies irrespective of zero cell counts.}
  \item{allstudies}{A logical indicating if studies with zero or all
    events in both groups are to be included in the meta-analysis
    (applies only if \code{sm} is equal to \code{"RR"} or \code{"OR"}).}
  \item{MH.exact}{A logical indicating if \code{incr} is not to be added
    to all cell frequencies for studies with a zero cell count to
    calculate the pooled estimate based on the Mantel-Haenszel method.}
  \item{RR.cochrane}{A logical indicating if 2*\code{incr} instead of
    1*\code{incr} is to be added to \code{n.e} and \code{n.c} in the
    calculation of the risk ratio (i.e., \code{sm="RR"}) for studies
    with a zero cell. This is used in RevMan 5, the
    Cochrane Collaboration's program for preparing and maintaining
    Cochrane reviews.}
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
    funnel plot asymmetry is to be used. Either \code{"rank"},
    \code{"linreg"}, \code{"mm"}, \code{"count"}, \code{"score"}, or
    \code{"peters"}, can be abbreviated. See function \code{\link{metabias}}}
  \item{backtransf}{A logical indicating whether results for odds
    ratio (\code{sm="OR"}) and risk ratio (\code{sm="RR"}) should be
    back transformed in printouts and plots. If TRUE (default),
    results will be presented as odds ratios and risk ratios;
    otherwise log odds ratios and log risk ratios will be shown.}
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
  \item{print.CMH}{A logical indicating whether result of the
    Cochran-Mantel-Haenszel test for overall effect should be printed.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in meta object.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if \code{incr} is added to studies with zero cell
    frequencies).}
}

\details{
  Treatment estimates and standard errors are calculated for each
  study. The following measures of treatment effect are available:
  \itemize{
   \item Risk ratio (\code{sm="RR"})
   \item Odds ratio (\code{sm="OR"})
   \item Risk difference (\code{sm="RD"})
   \item Arcsine difference (\code{sm="ASD"})
  }

  For studies with a zero cell count, by default, 0.5 is added to all
  cell frequencies of these studies; if \code{incr} is \code{"TACC"} a
  treatment arm continuity correction is used instead (Sweeting et
  al., 2004; Diamond et al., 2007). Treatment estimates and standard
  errors are only calculated for studies with zero or all events in
  both groups if \code{allstudies} is \code{TRUE}.
  
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

  By default, both fixed effect and random effects models are
  considered (see arguments \code{comb.fixed} and
  \code{comb.random}). If \code{method} is \code{"MH"} (default), the
  Mantel-Haenszel method is used to calculate the fixed effect
  estimate; if \code{method} is \code{"Inverse"}, inverse variance
  weighting is used for pooling; finally, if \code{method} is
  \code{"Peto"}, the Peto method is used for pooling. By default, the
  DerSimonian-Laird estimate (1986) is used in the random effects
  model (\code{method.tau="DL"}). For the Peto method, Peto's log odds
  ratio, i.e. \code{(O-E)/V} and its standard error \code{sqrt(1/V)}
  with \code{O-E} and \code{V} denoting "Observed minus Expected" and
  "V", are utilised in the random effects model. Accordingly, results
  of a random effects model using \code{sm="Peto"} can be (slightly)
  different to results from a random effects model using
  \code{sm="MH"} or \code{sm="Inverse"}.
  
  For the Mantel-Haenszel method, by default (if \code{MH.exact} is
  FALSE), 0.5 is added to all cell frequencies of a study with a zero cell
  count in the calculation of the pooled risk ratio or odds ratio as
  well as the estimation of the variance of the pooled risk difference,
  risk ratio or odds ratio. This approach is also used in other
  software, e.g. RevMan 5 and the Stata procedure metan. According to
  Fleiss (in Cooper & Hedges, 1994), there is no need to add 0.5 to a
  cell frequency of zero to calculate the Mantel-Haenszel estimate and
  he advocates the exact method (\code{MH.exact}=TRUE). Note, the
  estimate based on the exact method is not defined if the number of
  events is zero in all studies either in the experimental or control
  group.

  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metabin object by only specifying
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
  An object of class \code{c("metabin", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{event.e, n.e, event.c, n.c, studlab,}{}
  \item{sm, method, incr, allincr, addincr, }{}
  \item{allstudies, MH.exact, RR.cochrane, warn,}{}
  \item{level, level.comb, comb.fixed, comb.random,}{}
  \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{}
  \item{tau.common, title, complab, outclab,}{}
  \item{label.e, label.c, label.left, label.right,}{}
  \item{byvar, bylab, print.byvar}{As defined above.}
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
  \item{Q.CMH}{Cochran-Mantel-Haenszel test statistic for overall effect.}
  \item{incr.e, incr.c}{Increment added to cells in the experimental and
    control group, respectively.}
  \item{sparse}{Logical flag indicating if any study included in
    meta-analysis has any zero cell frequencies.}
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
  \item{n.e.w}{Number of observations in experimental group in
    subgroups - if \code{byvar} is not missing.}
  \item{event.c.w}{Number of events in control group in subgroups - if
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
  \item{call}{Function call.}
  \item{version}{Version of R package meta used to create object.}
}

\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.

  Diamond GA, Bax L, Kaul S (2007),
  Uncertain Effects of Rosiglitazone on the Risk for Myocardial
  Infarction and Cardiovascular Death.
  \emph{Annals of Internal Medicine}, \bold{147}, 578--581.

  DerSimonian R & Laird N (1986),
  Meta-analysis in clinical trials. \emph{Controlled Clinical Trials},
  \bold{7}, 177--188.

  Fleiss JL (1993),
  The statistical basis of meta-analysis.
  \emph{Statistical Methods in Medical Research}, \bold{2}, 121--145.

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
  
  \emph{Review Manager (RevMan)} [Computer program]. Version
  5.3. Copenhagen: The Nordic Cochrane Centre, The Cochrane
  Collaboration, 2014.
  
  Paule RC & Mandel J (1982),
  Consensus values and weighting factors.
  \emph{Journal of Research of the National Bureau of Standards},
  \bold{87}, 377--385.
  
  Pettigrew HM, Gart JJ, Thomas DG (1986),
  The bias and higher cumulants of the logarithm of a binomial
  variate. \emph{Biometrika}, \bold{73}, 425--435.

  Ruecker G, Schwarzer G, Carpenter JR (2008)
  Arcsine test for publication bias in meta-analyses with binary
  outcomes. \emph{Statistics in Medicine}, \bold{27}, 746--763.
    
  StataCorp. 2011.
  \emph{Stata Statistical Software: Release 12}. College Station, TX:
  StataCorp LP.
  
  Sweeting MJ, Sutton AJ, Lambert PC (2004),
  What to add to nothing? Use and avoidance of continuity corrections
  in meta-analysis of sparse data.
  \emph{Statistics in Medicine}, \bold{23}, 1351--1375.

  Viechtbauer W (2010),
  Conducting Meta-Analyses in R with the Metafor Package.
  \emph{Journal of Statistical Software}, \bold{36}, 1--48.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{update.meta}}, \code{\link{funnel}}, \code{\link{metabias}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
metabin(10, 20, 15, 20, sm="OR", warn=FALSE)

#
# Different results:
#
metabin(0, 10, 0, 10, sm="OR", warn=FALSE)
metabin(0, 10, 0, 10, sm="OR", allstudies=TRUE, warn=FALSE)


data(Olkin95)

meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", method="I")
summary(meta1)
funnel(meta1)

meta2 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=Olkin95$year<1970,
                 sm="RR", method="I")
summary(meta2)
forest(meta2)

meta3 <- update(meta2, sm="OR")
summary(meta3)

meta4 <- update(meta3, method="MH")
summary(meta4)

meta5 <- update(meta3, method="Peto")
summary(meta5)
}
\keyword{htest}
