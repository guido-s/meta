\name{metamean}

\alias{metamean}

\title{Meta-analysis of single means}

\description{
  Calculation of an overall mean from studies reporting a single mean
  using the inverse varinace method for pooling; inverse variance
  weighting is used for pooling.
}

\usage{
metamean(n, mean, sd, studlab,
         data=NULL, subset=NULL, exclude=NULL,
         sm=gs("smmean"),
         level=gs("level"), level.comb=gs("level.comb"),
         comb.fixed=gs("comb.fixed"), comb.random=gs("comb.random"),
         hakn=gs("hakn"),
         method.tau=gs("method.tau"), tau.preset=NULL, TE.tau=NULL,
         tau.common=gs("tau.common"),
         prediction=gs("prediction"), level.predict=gs("level.predict"),
         null.effect=NA,
         method.bias=gs("method.bias"),
         backtransf=gs("backtransf"),
         title=gs("title"), complab=gs("complab"), outclab="",
         byvar, bylab, print.byvar=gs("print.byvar"),
         byseparator=gs("byseparator"),
         keepdata=gs("keepdata"),
         warn=gs("warn"),
         control=NULL)
}

\arguments{
  \item{n}{Number of observations.}
  \item{mean}{Estimated mean.}
  \item{sd}{Standard deviation.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be
    used.}
  \item{exclude}{An optional vector specifying studies to exclude from
    meta-analysis, however, to include in printouts and forest plots.}
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
  \item{hakn}{A logical indicating whether the method by Hartung and
    Knapp should be used to adjust test statistics and confidence
    intervals.}	  
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
    back transformed in printouts and plots for \code{sm="MLN"}. If
    TRUE (default), results will be presented as means; otherwise
    logarithm of means will be shown.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{sm}{A character string indicating which summary measure
    (\code{"MRAW"} or \code{"MLN"}) is to be used for pooling of
    studies.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{n}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{byseparator}{A character string defining the separator between
    label and levels of grouping variable.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in meta object.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if studies are excluded from meta-analysis due to zero
    standard deviations).}
  \item{control}{An optional list to control the iterative process to
    estimate the between-study variance tau^2. This argument is passed
    on to \code{\link[metafor]{rma.uni}}.}  
}

\details{
  Fixed effect and random effects meta-analysis of single means to
  calculate an overall mean; inverse variance weighting is used for
  pooling. The following transformations of means are implemented to
  calculate an overall mean:
  \itemize{
  \item Raw, i.e. untransformed, means (\code{sm="MRAW"}, default)
  \item Log transformed means (\code{sm="MLN"})
  }
     
  Note, you should use R function \code{\link{metacont}} to compare
  means of pairwise comparisons instead of using \code{metamean} for
  each treatment arm separately which will break randomisation in
  randomised controlled trials.
  
  Calculations are conducted on the log scale if
  \code{sm="ROM"}. Accordingly, list elements \code{TE},
  \code{TE.fixed}, and \code{TE.random} contain the logarithm of
  means. In printouts and plots these values are back transformed if
  argument \code{backtransf=TRUE}.
  
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

  The function \code{metagen} is called internally to calculate
  individual and overall treatment estimates and standard errors.

  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metamean object by only specifying
  arguments which should be changed.

  For the random effects, the method by Hartung and Knapp (2001) /
  Knapp and Hartung (2003) is used to adjust test statistics and
  confidence intervals if argument \code{hakn=TRUE}.

  The DerSimonian-Laird estimate (1986) is used in the random effects
  model if \code{method.tau="DL"}. The iterative Paule-Mandel method
  (1982) to estimate the between-study variance is used if argument
  \code{method.tau="PM"}. Internally, R function \code{paulemandel} is
  called which is based on R function \code{mpaule.default} from R
  package \bold{metRology} from S.L.R. Ellison <s.ellison at
  lgc.co.uk>.
  
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
  An object of class \code{c("metamean", "meta")} with corresponding
  \code{print}, \code{summary}, and \code{forest} functions. The
  object is a list containing the following components:
  \item{n, mean, sd,}{}
  \item{studlab, exclude, sm, level, level.comb,}{}
  \item{comb.fixed, comb.random,}{}
  \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{}
  \item{tau.common, title, complab, outclab,}{}
  \item{byvar, bylab, print.byvar, byseparator, warn}{As defined above.}
  \item{TE, seTE}{Estimated effect (mean or log mean) and standard
    error of individual studies.}
  \item{lower, upper}{Lower and upper confidence interval limits
    for individual studies.}
  \item{zval, pval}{z-value and p-value for test of overall effect for
    individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall effect (mean or log
    mean) and standard error (fixed effect model).}  
  \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
    limits (fixed effect model).}
  \item{zval.fixed, pval.fixed}{z-value and p-value for test of
    overall effect (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall effect (mean or log
    mean) and standard error (random effects model).}  
  \item{lower.random, upper.random}{Lower and upper confidence
    interval limits (random effects model).}
  \item{zval.random, pval.random}{z-value or t-value and corresponding
    p-value for test of overall effect (random effects model).}
  \item{prediction, level.predict}{As defined above.}
  \item{seTE.predict}{Standard error utilised for prediction interval.}
  \item{lower.predict, upper.predict}{Lower and upper limits of
    prediction interval.}
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
  \item{TE.fixed.w, seTE.fixed.w}{Estimated effect and standard error
    in subgroups (fixed effect model) - if \code{byvar} is not
    missing.}
  \item{lower.fixed.w, upper.fixed.w}{Lower and upper confidence
    interval limits in subgroups (fixed effect model) - if
    \code{byvar} is not missing.}
  \item{zval.fixed.w, pval.fixed.w}{z-value and p-value for test of
    treatment effect in subgroups (fixed effect model) - if
    \code{byvar} is not missing.}
  \item{TE.random.w, seTE.random.w}{Estimated effect and standard
    error in subgroups (random effects model) - if \code{byvar} is not
    missing.}
  \item{lower.random.w, upper.random.w}{Lower and upper confidence
    interval limits in subgroups (random effects model) - if
    \code{byvar} is not missing.}
  \item{zval.random.w, pval.random.w}{z-value or t-value and
    corresponding p-value for test of effect in subgroups (random
    effects model) - if \code{byvar} is not missing.}
  \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
    random effects model) - if \code{byvar} is not missing.}
  \item{df.hakn.w}{Degrees of freedom for test of effect for
    Hartung-Knapp method in subgroups - if \code{byvar} is not missing
    and \code{hakn=TRUE}.}
  \item{n.e.w}{Number of observations in experimental group in
    subgroups - if \code{byvar} is not missing.}
  \item{n.c.w}{Number of observations in control group in subgroups -
    if \code{byvar} is not missing.}
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
  \item{keepdata}{As defined above.}
  \item{data}{Original data (set) used in function call (if
    \code{keepdata=TRUE}).}
  \item{subset}{Information on subset of original data used in
    meta-analysis (if \code{keepdata=TRUE}).}
  \item{call}{Function call.}
  \item{version}{Version of R package \bold{meta} used to create object.}
}

\references{
  DerSimonian R & Laird N (1986),
  Meta-analysis in clinical trials. \emph{Controlled Clinical Trials},
  \bold{7}, 177--88.

  Hartung J & Knapp G (2001),
  On tests of the overall treatment effect in meta-analysis with
  normally distributed responses.
  \emph{Statistics in Medicine}, \bold{20}, 1771--82.
  doi: 10.1002/sim.791 .

  Higgins JPT, Thompson SG, Spiegelhalter DJ (2009),
  A re-evaluation of random-effects meta-analysis.
  \emph{Journal of the Royal Statistical Society: Series A},
  \bold{172}, 137--59.

  Knapp G & Hartung J (2003),
  Improved Tests for a Random Effects Meta-regression with a Single
  Covariate.
  \emph{Statistics in Medicine}, \bold{22}, 2693--710,
  doi: 10.1002/sim.1482 .
  
  Paule RC & Mandel J (1982),
  Consensus values and weighting factors.
  \emph{Journal of Research of the National Bureau of Standards},
  \bold{87}, 377--85.

  Viechtbauer W (2010),
  Conducting Meta-Analyses in R with the Metafor Package.
  \emph{Journal of Statistical Software}, \bold{36}, 1--48.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{update.meta}}, \code{\link{metamean}}, \code{\link{metagen}}}

\examples{
m1 <- metamean(rep(100, 3), 1:3, rep(1, 3))
m2 <- update(m1, sm="MLN")

m1
m2

# With test for overall mean equal to 2
#
update(m1, null.effect=2)
update(m2, null.effect=2)

# Print results without back-transformation
#
print(m1, backtransf=FALSE)
update(m2, backtransf=FALSE)
update(m1, null.effect=2, backtransf=FALSE)
update(m2, null.effect=2, backtransf=FALSE)
}

\keyword{htest}
