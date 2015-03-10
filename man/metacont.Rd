\name{metacont}

\alias{metacont}

\title{Meta-analysis of continuous outcome data}

\description{
  Calculation of fixed and random effects estimates for meta-analyses
  with continuous outcome data; inverse variance weighting is used for
  pooling.
}

\usage{
metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab,
         data=NULL, subset=NULL,
         sm=.settings$smcont, pooledvar=.settings$pooledvar,
	 method.smd=.settings$method.smd, sd.glass=.settings$sd.glass,
	 exact.smd=.settings$exact.smd,
         level=.settings$level, level.comb=.settings$level.comb,
         comb.fixed=.settings$comb.fixed, comb.random=.settings$comb.random,
         hakn=.settings$hakn,
         method.tau=.settings$method.tau, tau.preset=NULL, TE.tau=NULL,
         tau.common=.settings$tau.common,
         prediction=.settings$prediction, level.predict=.settings$level.predict,
         method.bias=.settings$method.bias,
         title=.settings$title, complab=.settings$complab, outclab="",
         label.e=.settings$label.e, label.c=.settings$label.c,
         label.left=.settings$label.left, label.right=.settings$label.right,
         byvar, bylab, print.byvar=.settings$print.byvar,
         keepdata=.settings$keepdata,
         warn=.settings$warn)
}

\arguments{
  \item{n.e}{Number of observations in experimental group.}
  \item{mean.e}{Estimated mean in experimental group.}
  \item{sd.e}{Standard deviation in experimental group.}
  \item{n.c}{Number of observations in control group.}
  \item{mean.c}{Estimated mean in control group.}
  \item{sd.c}{Standard deviation in control group.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be
  used.}
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
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{label.e}{Label for experimental group.}
  \item{label.c}{Label for control group.}
  \item{label.left}{Graph label on left side of forest plot.}
  \item{label.right}{Graph label on right side of forest plot.}
  \item{sm}{A character string indicating which summary measure
    (\code{"MD"} or \code{"SMD"}) is to be used for pooling of
    studies.}
  \item{pooledvar}{A logical indicating if a pooled variance should be
    used for the mean difference (\code{sm="MD"}).}
  \item{method.smd}{A character string indicating which method is used
    to estimate the standardised mean difference
    (\code{sm="SMD"}). Either \code{"Hedges"} for Hedges' g (default),
    \code{"Cohen"} for Cohen's d, or \code{"Glass"} for Glass' delta,
    can be abbreviated.}
  \item{sd.glass}{A character string indicating which standard
    deviation is used in the denominator for Glass' method to estimate
    the standardised mean difference. Either \code{"control"} using
    the standard deviation in the control group (\code{sd.c}) or
    \code{"experimental"} using the standard deviation in the
    experimental group (\code{sd.e}), can be abbreviated.}
  \item{exact.smd}{A logical indicating whether exact formulae should
    be used in estimation of the standardised mean difference and its
    standard error (see Details).}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{n.e}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in meta object.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if studies are excluded from meta-analysis due to zero
    standard deviations).}
}

\details{
  Calculation of fixed and random effects estimates for meta-analyses
  with continuous outcome data; inverse variance weighting is used for
  pooling. By default, the DerSimonian-Laird estimate (1986) is used
  in the random effects model (\code{method.tau="DL"}).

  Two different types of summary measures are available for continuous
  outcomes: the mean difference (argument \code{sm="MD"}) and the
  standardised mean difference (\code{sm="SMD"}). For the standardised
  mean difference three methods are implemented:
  \itemize{
   \item Hedges' g (default, \code{method.smd="Hedges"}) - see Hedges (1981)
   \item Cohen's d (\code{method.smd="Cohen"}) - see Cohen (1988)
   \item Glass' delta (\code{method.smd="Glass"}) - see Glass (1976)
  }
  Hedges (1981) calculated the exact bias in Cohen's d which is a
  ratio of gamma distributions with the degrees of freedom, i.e. total
  sample size minus two, as argument. By default (argument
  \code{exact.smd=FALSE}), an accurate approximation of this bias
  provided in Hedges (1981) is utilised for Hedges' g as well as its
  standard error; these approximations are also used in RevMan
  5. Following Borenstein et al. (2009) these approximations are not
  used in the estimation of Cohen's d. White and Thomas (2005) argued
  that approximations are unnecessary with modern software and
  accordingly promote to use the exact formulae; this is possible
  using argument \code{exact.smd=TRUE}. For Hedges' g the exact
  formulae are used to calculate the standardised mean difference as
  well as the standard error; for Cohen's d the exact formula is only
  used to calculate the standard error. In typical applications (with
  sample sizes above 10), the differences between using the exact
  formulae and the approximation will be minimal.
  
  For Glass' delta, by default (argument \code{sd.glass="control"}),
  the standard deviation in the control group (\code{sd.c}) is used in
  the denominator of the standard mean difference. The standard
  deviation in the experimental group (\code{sd.e}) can be used by
  specifying \code{sd.glass="experimental"}.
  
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

  The function \code{metagen} is called internally to calculate
  individual and overall treatment estimates and standard errors.

  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metacont object by only specifying
  arguments which should be changed.

  For the random effects, the method by Hartung and Knapp (2003) is
  used to adjust test statistics and confidence intervals if argument
  \code{hakn=TRUE}.

  The iterative Paule-Mandel method (1982) to estimate the
  between-study variance is used if argument
  \code{method.tau="PM"}. Internally, R function \code{paulemandel} is
  called which is based on R function \code{mpaule.default} from R
  package metRology from S.L.R. Ellison <s.ellison at lgc.co.uk>.
  
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
  An object of class \code{c("metacont", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{n.e, mean.e, sd.e,}{}
  \item{n.c, mean.c, sd.c,}{}
  \item{studlab, sm, level, level.comb,}{}
  \item{comb.fixed, comb.random,}{}
  \item{pooledvar, method.smd, sd.glass,}{}
  \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{}
  \item{tau.common, title, complab, outclab,}{}
  \item{label.e, label.c, label.left, label.right,}{}
  \item{byvar, bylab, print.byvar, warn}{As defined above.}
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
  \item{Q}{Heterogeneity statistic.}
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
  \item{version}{Version of R package meta used to create object.}
}

\references{
  Borenstein et al. (2009),
  \emph{Introduction to Meta-Analysis},
  Chichester: Wiley.

  Cohen J (1988),
  Statistical Power Analysis for the Behavioral Sciences (second ed.).
  Lawrence Erlbaum Associates.

  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.

  DerSimonian R & Laird N (1986),
  Meta-analysis in clinical trials. \emph{Controlled Clinical Trials},
  \bold{7}, 177--188.

  Glass G (1976),
  Primary, secondary, and meta-analysis of research.
  \emph{Educational Researcher}, \bold{5}, 3--8.

  Hartung J & Knapp G (2001),
  On tests of the overall treatment effect in meta-analysis with
  normally distributed responses.
  \emph{Statistics in Medicine}, \bold{20}, 1771--82.
  doi: 10.1002/sim.791 .
  
  Hedges LV (1981),
  Distribution theory for Glass's estimator of effect size and related
  estimators.
  \emph{Journal of Educational and Behavioral Statistics},
  \bold{6}, 107--28.

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
  
  \emph{Review Manager (RevMan)} [Computer program]. Version
  5.3. Copenhagen: The Nordic Cochrane Centre, The Cochrane
  Collaboration, 2014.

  Viechtbauer W (2010),
  Conducting Meta-Analyses in R with the Metafor Package.
  \emph{Journal of Statistical Software}, \bold{36}, 1--48.

  White IR, Thomas J (2005),
  Standardized mean differences in individually-randomized and
  cluster-randomized trials, with applications to meta-analysis.
  \emph{Clinical Trials}, \bold{2}, 141--51.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{update.meta}}, \code{\link{metabin}}, \code{\link{metagen}}}

\examples{
data(Fleiss93cont)
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont, sm="SMD")
meta1
forest(meta1)

meta2 <- metacont(Fleiss93cont$n.e, Fleiss93cont$mean.e,
                  Fleiss93cont$sd.e,
                  Fleiss93cont$n.c, Fleiss93cont$mean.c,
                  Fleiss93cont$sd.c,
                  sm="SMD")
meta2

data(amlodipine)
meta3 <- metacont(n.amlo, mean.amlo, sqrt(var.amlo),
                  n.plac, mean.plac, sqrt(var.plac),
                  data=amlodipine, studlab=study)
summary(meta3)

# Use pooled variance
#
meta4 <- metacont(n.amlo, mean.amlo, sqrt(var.amlo),
                  n.plac, mean.plac, sqrt(var.plac),
                  data=amlodipine, studlab=study,
                  pooledvar=TRUE)
summary(meta4)

# Use Cohen's d instead of Hedges' g as effect measure
#
meta5 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont,
                  sm="SMD", method.smd="Cohen")
meta5

# Use Glass' delta instead of Hedges' g as effect measure
#
meta6 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont,
                  sm="SMD", method.smd="Glass")
meta6

# Use Glass' delta based on the standard deviation in the experimental group
#
meta7 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont,
                  sm="SMD", method.smd="Glass", sd.glass="experimental")
meta7

# Calculate Hedges' g based on exact formulae
#
update(meta1, exact.smd=TRUE)
}

\keyword{htest}
