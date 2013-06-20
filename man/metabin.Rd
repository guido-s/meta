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
        data=NULL, subset=NULL, method="MH",
        sm=ifelse(!is.na(charmatch(method, c("Peto", "peto"),
                                   nomatch=NA)), "OR", "RR"),
        incr=0.5, allincr=FALSE, addincr=FALSE, allstudies=FALSE,
        MH.exact=FALSE, RR.cochrane=FALSE,
        level=0.95, level.comb=level,
        comb.fixed=TRUE, comb.random=TRUE,
        hakn=FALSE,
        method.tau="DL", tau.preset=NULL, TE.tau=NULL,
        tau.common=FALSE,
        prediction=FALSE, level.predict=level,
        method.bias=NULL,
        title="", complab="", outclab="",
        label.e="Experimental", label.c="Control",
        label.left="", label.right="",
        byvar, bylab, print.byvar=TRUE,
        print.CMH=FALSE, warn=TRUE)
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
    (\code{"RR"}, \code{"OR"}, \code{"RD"}, or \code{"AS"}) is to be used
    for pooling of studies, see Details.}
  \item{incr}{Could be either a numerical value which is added to each cell
  frequency for studies with a zero cell count or the character string
  \code{"TA"} which stands for treatment arm continuity correction, see
  Details.}
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
    \code{"DL"}, \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
    \code{"HE"}, or \code{"EB"}, can be abbreviated.}
  \item{tau.preset}{Prespecified value for between-study variance \eqn{\tau^2}.}
  \item{TE.tau}{Overall treatment effect used to estimate the
  between-study variance \eqn{\tau^2}.}
  \item{tau.common}{A logical indicating whether tau-squared should be
    the same across subgroups.}
  \item{method.bias}{A character string indicating which test for
    funnel plot asymmetry is to be used. Either \code{"rank"},
    \code{"linreg"}, \code{"mm"}, \code{"count"}, \code{"score"}, or
    \code{"peters"}, can be abbreviated.}
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
   \item Arcsine difference (\code{sm="AS"})
  }

  For studies with a zero cell count, by default, 0.5 is added to
  all cell frequencies of these studies; if \code{incr} is
  \code{"TA"} a treatment arm continuity correction is used instead
  (Sweeting et al., 2004; Diamond et al., 2007). Treatment estimates
  and standard errors are only calculated for studies with zero or all
  events in both groups if \code{allstudies} is \code{TRUE}.
  
  Internally, both fixed effect and random effects models are calculated
  regardless of values choosen for arguments \code{comb.fixed} and
  \code{comb.random}. Accordingly, the estimate for the random effects
  model can be extracted from component \code{TE.random} of an object
  of class \code{"meta"} even if \code{comb.random=FALSE}. However, all
  functions in R package \code{meta} will adequately consider the values
  for \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.meta}} will not print results for the random effects
  model if \code{comb.random=FALSE}.

  By default, both fixed effect and random effects models are
  considered (arguments \code{comb.fixed=TRUE} and
  \code{comb.random=TRUE}). If \code{method} is \code{"MH"} (default),
  the Mantel-Haenszel method is used to calculate the fixed effect
  estimate; if \code{method} is \code{"Inverse"}, inverse variance
  weighting is used for pooling; finally, if \code{method} is
  \code{"Peto"}, the Peto method is used for pooling. By default, the
  DerSimonian-Laird estimate is used in the random effects model (see
  paragraph on argument \code{method.tau}). For the Peto method,
  Peto's log odds ratio, i.e. \code{(O-E)/V} and its standard error
  \code{sqrt(1/V)} with \code{O-E} and \code{V} denoting "Observed
  minus Expected" and "V", are utilised in the random effects
  model. Accordingly, results of a random effects model using
  \code{sm="Peto"} can be (slightly) different to results from a
  random effects model using \code{sm="MH"} or \code{sm="Inverse"}.
  
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

  If R package metafor (Viechtbauer 2010) is installed, the following
  statistical methods are also available.

  For the random effects model (argument \code{comb.random=TRUE}), the
  method by Hartung and Knapp (Hartung, Knapp 2001; Knapp, Hartung
  2003) is used to adjust test statistics and confidence intervals if
  argument \code{hakn=TRUE} (internally R function \code{rma.uni} of R
  package metafor is called).

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
  \item{tau}{Square-root of between-study variance.}
  \item{se.tau}{Standard error of square-root of between-study variance.}
  \item{C}{Scaling factor utilised internally to calculate common
    tau-squared across subgroups.}
  \item{Q.CMH}{Cochran-Mantel-Haenszel test statistic for overall effect.}
  \item{incr.e, incr.c}{Increment added to cells in the experimental and
    control group, respectively}
  \item{sparse}{Logical flag indicating if any study included in
    meta-analysis has any zero cell frequencies.}
  \item{df.hakn}{Degrees of freedom for test of treatment effect for
    Hartung-Knapp method (only if \code{hakn=TRUE}).}
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
  \bold{172}, 137-159.

  Knapp G & Hartung J (2003),
  Improved Tests for a Random Effects Meta-regression with a Single
  Covariate.
  \emph{Statistics in Medicine}, \bold{22}, 2693-710,
  doi: 10.1002/sim.1482 .
  
  \emph{Review Manager (RevMan)} [Computer program]. Version
  5.1. Copenhagen: The Nordic Cochrane Centre, The Cochrane
  Collaboration, 2011.

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

\seealso{\code{\link{funnel}}, \code{\link{metabias}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
metabin(10, 20, 15, 20, sm="OR", warn=FALSE)

##
## Different results:
##
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
}
\keyword{htest}
