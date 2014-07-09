\name{metaprop}
\alias{metaprop}

\title{Meta-analysis of single proportions}

\description{
  Calculation of an overall proportion from studies reporting a single
  proportion.
}

\usage{
metaprop(event, n, studlab,
         data=NULL, subset=NULL,
         sm=.settings$smprop,
         incr=.settings$incr, allincr=.settings$allincr,
         addincr=.settings$addincr,
         method.ci=.settings$method.ci,
         level=.settings$level, level.comb=.settings$level.comb,
         comb.fixed=.settings$comb.fixed, comb.random=.settings$comb.random,
         hakn=.settings$hakn,
         method.tau=.settings$method.tau, tau.preset=NULL, TE.tau=NULL,
         tau.common=.settings$tau.common,
         prediction=.settings$prediction, level.predict=.settings$level.predict,
         method.bias=.settings$method.bias,
         title=.settings$title, complab=.settings$complab, outclab="",
         byvar, bylab, print.byvar=.settings$print.byvar,
         keepdata=.settings$keepdata,
         warn=.settings$warn)
}

\arguments{
  \item{event}{Number of events.}
  \item{n}{Number of observations.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., event and n.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating which summary measure
    (\code{"PFT"}, \code{"PAS"}, \code{"PRAW"}, \code{"PLN"}, or
    \code{"PLOGIT"}) is to be used for pooling of studies, see Details.}
  \item{incr}{A numeric which is added to each cell frequency for
    studies with a zero cell count.}
  \item{allincr}{A logical indicating if \code{incr} is added to each
    cell frequency of all studies if at least one study has a zero cell
    count. If FALSE (default), \code{incr} is added only to each cell frequency of
    studies with a zero cell count.}
  \item{addincr}{A logical indicating if \code{incr} is added to each cell
    frequency of all studies irrespective of zero cell counts.}
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
  \item{warn}{A logical indicating whether the addition of \code{incr}
    to studies with zero cell frequencies should result in a warning.}
}

\details{
  Fixed effect and random effects meta-analysis of single proportions
  to calculate an overall proportion. By default, the
  DerSimonian-Laird estimate (1986) is used in the random effects
  model (\code{method.tau="DL"}).
  
  The following transformations of proportions are implemented to
  calculate an overall proportion:

  \itemize{

  \item Logit transformation (\code{sm="PLOGIT"}, default)
  
  \item Log transformation (\code{sm="PLN"})

  \item Freeman-Tukey Double arcsine transformation (\code{sm="PFT"})  
  
  \item Arcsine transformation (\code{sm="PAS"})
  
  \item Raw, i.e. untransformed, proportions (\code{sm="PRAW"})}

  If the summary measure is equal to "PRAW", "PLN", or "PLOGIT", a
  continuity correction is applied if any studies has a zero cell
  count. By default, 0.5 is added to all cell frequencies of studies
  with a zero cell count (argument \code{incr}).

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
  calculation disregarding the chosen summary measure.
  
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
  meta-analysis of an existing metaprop object by only specifying
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
  An object of class \code{c("metaprop", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{event, n, studlab,}{}
  \item{sm, incr, allincr, addincr, method.ci,}{}
  \item{level, level.comb,}{As defined above.}
  \item{comb.fixed, comb.random,}{}
  \item{hakn, method.tau, tau.preset, TE.tau, method.bias,}{}
  \item{tau.common, title, complab, outclab,}{}
  \item{byvar, bylab, print.byvar, warn}{}
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
  \item{sm}{A character string: \code{"proportion"}}
  \item{method}{A character string indicating method used
    for pooling: \code{"Inverse"}}
  \item{df.hakn}{Degrees of freedom for test of treatment effect for
  Hartung-Knapp method (only if \code{hakn=TRUE}).}
  \item{incr.event}{Increment added to number of events.}
  \item{keepdata}{As defined above.}
  \item{data}{Original data (set) used in function call (if
    \code{keepdata=TRUE}).}
  \item{subset}{Information on subset of original data used in
    meta-analysis (if \code{keepdata=TRUE}).}
  \item{call}{Function call.}
  \item{version}{Version of R package meta used to create object.}
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

  Viechtbauer W (2010),
  Conducting Meta-Analyses in R with the Metafor Package.
  \emph{Journal of Statistical Software}, \bold{36}, 1--48.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{update.meta}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
metaprop(4:1, c(10, 20, 30, 40))
metaprop(4:1, c(10, 20, 30, 40), sm="PAS")
metaprop(4:1, c(10, 20, 30, 40), sm="PRAW")
metaprop(4:1, c(10, 20, 30, 40), sm="PLN")
metaprop(4:1, c(10, 20, 30, 40), sm="PFT")

forest(metaprop(4:1, c(10, 20, 30, 40)))


m1 <- metaprop(c(0, 0, 10, 10), rep(100, 4))
m2 <- metaprop(c(0, 0, 10, 10), rep(100, 4), incr=0.1)

summary(m1)
summary(m2)

forest(m1)
forest(m2)

## Example from Miller (1978):
death <- c(3, 6, 10, 1)
animals <- c(11, 17, 21, 6)
##
m3 <- metaprop(death, animals, sm="PFT")
forest(m3)

## Data examples from Newcombe (1998)
##
event <- c(81, 15, 0, 1)
n <- c(263, 148, 20, 29)
##
m1 <- metaprop(event, n, sm="PLOGIT", method.ci="SA")
m2 <- update(m1, method.ci="SACC")
m3 <- update(m1, method.ci="WS")
m4 <- update(m1, method.ci="WSCC")
m5 <- update(m1, method.ci="CP")
##
lower <- round(rbind(NA, m1$lower, m2$lower, NA, m3$lower, m4$lower, NA, m5$lower), 4)
upper <- round(rbind(NA, m1$upper, m2$upper, NA, m3$upper, m4$upper, NA, m5$upper), 4)
##
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
##
## Newcombe (1998), Table I, methods 1-5:
##
tab1

## Same confidence interval, i.e. unaffected by choice of summary measure
##
print(metaprop(event, n, sm="PLOGIT", method.ci="WS"), ma=FALSE)
print(metaprop(event, n, sm="PLN", method.ci="WS"), ma=FALSE)
print(metaprop(event, n, sm="PFT", method.ci="WS"), ma=FALSE)
print(metaprop(event, n, sm="PAS", method.ci="WS"), ma=FALSE)
print(metaprop(event, n, sm="PRAW", method.ci="WS"), ma=FALSE)
}

\keyword{htest}
