\name{metainc}

\alias{metainc}

\title{Meta-analysis of incidence rates}

\description{
  Calculation of fixed and random effects estimates (incidence rate
  ratio or incidence rate difference) for meta-analyses with event
  counts. Mantel-Haenszel, Cochran, and inverse variance method are
  available for pooling.
}

\usage{
metainc(event.e, time.e, event.c, time.c, studlab,
        data=NULL, subset=NULL, method="MH",
        sm="IRR",
        incr=0.5, allincr=FALSE, addincr=FALSE,
        level=0.95, level.comb=level,
        comb.fixed=TRUE, comb.random=TRUE,
        hakn=FALSE,
        method.tau="DL", tau.preset=NULL, TE.tau=NULL,
        tau.common=FALSE,
        prediction=FALSE, level.predict=level,
        method.bias="linreg",
        n.e=NULL, n.c=NULL,
        title="", complab="", outclab="",
        label.e="Experimental", label.c="Control",
        label.left="", label.right="",
        byvar, bylab, print.byvar=TRUE,
        keepdata=TRUE,
        warn=TRUE)
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
  \item{method}{A character string indicating which method is to be used
    for pooling of studies. One of \code{"MH"}, \code{"Inverse"}, or
    \code{"Cochran"}, can be abbreviated.}
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
    funnel plot asymmetry is to be used. Either \code{"linreg"} or
    \code{"rank"}, can be abbreviated.}
  \item{n.e}{Number of observations in experimental group (optional).}
  \item{n.c}{Number of observations in control group (optional).}
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
   \item Incidence Rate Ratio (\code{sm="IRR"})
   \item Incidence Rate Difference (\code{sm="IRD"})
  }

  For studies with a zero cell count, by default, 0.5 is added to
  all cell frequencies of these studies (argument \code{incr}).
  
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
  estimate (Greenland & Robbins, 1985); if \code{method} is
  \code{"Inverse"}, inverse variance weighting is used for pooling;
  finally, if \code{method} is \code{"Cochran"}, the Cochran method is
  used for pooling (Bayne-Jones, 1964, Chapter 8). By default, the
  DerSimonian-Laird estimate is used in the random effects model (see
  paragraph on argument \code{method.tau}).
  
  For Mantel-Haenszel and Cochran method, nothing is added to zero
  cell counts. Accordingly, Mantel-Haenszel and Cochran estimate are
  not defined if the number of events is zero in all studies either in
  the experimental or control group.

  A prediction interval for treatment effect of a new study is
  calculated (Higgins et al., 2009) if arguments \code{prediction} and
  \code{comb.random} are \code{TRUE}.

  R function \code{\link{update.meta}} can be used to redo the
  meta-analysis of an existing metainc object by only specifying
  arguments which should be changed.

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
}

\value{
  An object of class \code{c("metainc", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{event.e, time.e, event.c, time.c, studlab,}{}
  \item{sm, method, incr, allincr, addincr, warn,}{}
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
  \item{keepdata}{As defined above.}
  \item{data}{Original data (set) used in function call (if
    \code{keepdata=TRUE}).}
  \item{subset}{Information on subset of original data used in
    meta-analysis (if \code{keepdata=TRUE}).}
  \item{call}{Function call.}
  \item{version}{Version of R package meta used to create object.}
}

\references{
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
  \bold{172}, 137-159.

  Knapp G & Hartung J (2003),
  Improved Tests for a Random Effects Meta-regression with a Single
  Covariate.
  \emph{Statistics in Medicine}, \bold{22}, 2693-710,
  doi: 10.1002/sim.1482 .

  Bayne-Jones S et al. (1964),
  Smoking and Health: Report of the Advisory Committee to the Surgeon
  General of the United States. U-23 Department of Health, Education,
  and Welfare. Public Health Service Publication No. 1103.
  \url{http://profiles.nlm.nih.gov/ps/retrieve/ResourceMetadata/NNBBMQ}

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
}
\keyword{htest}
