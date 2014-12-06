\name{trimfill}

\alias{trimfill}
\alias{trimfill.meta}
\alias{trimfill.default}

\title{Trim-and-fill method to adjust for bias in meta-analysis}

\description{
  Trim-and-fill method for estimating and adjusting for the number and
  outcomes of missing studies in a meta-analysis.
}

\usage{
trimfill(x, ...)

\method{trimfill}{default}(x, seTE, left=NULL, ma.fixed=TRUE, type="L", n.iter.max=50,
         sm=NULL, studlab=NULL, level=0.95, level.comb=level,
         comb.fixed=FALSE, comb.random=TRUE,
         hakn=FALSE, method.tau="DL",
         prediction=FALSE, level.predict=level,
         backtransf=TRUE,
         silent=TRUE, ...)

\method{trimfill}{meta}(x, left=NULL, ma.fixed=TRUE, type="L", n.iter.max=50,
         level=x$level, level.comb=x$level.comb,
         comb.fixed=FALSE, comb.random=TRUE,
         hakn=x$hakn, method.tau=x$method.tau,
         prediction=x$prediction, level.predict=x$level.predict,
         backtransf=x$backtransf,
         silent=TRUE, ...)
}

\arguments{
  \item{x}{An object of class \code{meta}, or estimated treatment
    effect in individual studies.}
  \item{seTE}{Standard error of estimated treatment effect.}
  \item{left}{A logical indicating whether studies are supposed to be
    missing on the left or right side of the funnel plot. If NULL, the
    linear regression test for funnel plot symmetry (i.e., function
    \code{metabias(..., method="linreg")}) is used to determine whether
    studies are missing on the left or right side.}
  \item{ma.fixed}{A logical indicating whether a fixed effect or random
  effects model is used to estimate the number of missing studies.}
  \item{type}{A character indicating which method is used to estimate
    the number of missing studies. Either \code{"L"} or \code{"R"}.}
  \item{n.iter.max}{Maximum number of iterations to estimate number of
    missing studies.}
  \item{sm}{An optional character string indicating underlying summary
    measure, e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
    \code{"MD"}, \code{"SMD"}; ignored if \code{x} is of class
    \code{meta}.}
  \item{studlab}{An optional vector with study labels; ignored if
    \code{x} is of class \code{meta}.}
  \item{level}{The level used to calculate confidence intervals for
    individual studies. If existing, \code{x$level} is used as value for
    \code{level}; otherwise 0.95 is used.}
  \item{level.comb}{The level used to calculate confidence interval for
    the pooled estimate. If existing, \code{x$level.comb} is used as
    value for \code{level.comb}; otherwise 0.95 is used.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{hakn}{A logical indicating whether the method by Hartung and
          Knapp should be used to adjust test statistics and
          confidence intervals.}
  \item{method.tau}{A character string indicating which method is used
    to estimate the between-study variance \eqn{\tau^2}. Either
    \code{"DL"}, \code{"PM"}, \code{"REML"}, \code{"ML"}, \code{"HS"},
    \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be abbreviated.}
  \item{prediction}{A logical indicating whether a prediction interval
    should be printed.}
  \item{level.predict}{The level used to calculate prediction interval
    for a new study.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and plots. If
    \code{backtransf=TRUE}, results for \code{sm="OR"} are printed as
    odds ratios rather than log odds ratios and results for
    \code{sm="ZCOR"} are printed as correlations rather than Fisher's
    z transformed correlations, for example.}
  \item{silent}{A logical indicating whether basic information on
    iterations shown.}
  \item{\dots}{other arguments}
}

\details{
  The trim-and-fill method (Duval, Tweedie 2000a, 2000b) can be used
  for estimating and adjusting for the number and outcomes of missing
  studies in a meta-analysis. The method relies on scrutiny of one
  side of a funnel plot for asymmetry assumed due to publication bias.

  Three different methods have been proposed originally to estimate
  the number of missing studies. Two of these methods (L- and
  R-estimator) have been shown to perform better in simulations, and
  are available in this R function (argument \code{type}).

  A fixed effect or random effects model can be used to estimate the
  number of missing studies (argument \code{ma.fixed}). Furthermore, a
  fixed effect and/or random effects model can be used to summaries
  study results (arguments \code{comb.fixed} and
  \code{comb.random}). Simulation results (Peters et al. 2007)
  indicate that the fixed-random model, i.e. using a fixed effect
  model to estimate the number of missing studies and a random effects
  model to summaries results, (i) performs better than the fixed-fixed
  model, and (ii) performs no worse than and marginally better in
  certain situations than the random-random model. Accordingly, the
  fixed-random model is the default.

  An empirical comparison of the trim-and-fill method and the Copas
  selection model (Schwarzer et al. 2010) indicates that the
  trim-and-fill method leads to excessively conservative inference in
  practice. The Copas selection model is available in R package copas.

  The function \code{metagen} is called internally.
}

\value{
  An object of class \code{c("metagen", "meta", "trimfill")}. The object is a
  list containing the following components:
  \item{studlab, sm, left, ma.fixed, type,}{}
  \item{n.iter.max, level, level.comb, level.predict,}{As defined above.}
  \item{comb.fixed, comb.random, prediction, hakn, method.tau}{}
  \item{TE, seTE}{Estimated treatment effect and standard error of individual studies.}
  \item{lower, upper}{Lower and upper confidence interval limits
    for individual studies.}
  \item{zval, pval}{z-value and p-value for test of treatment
    effect for individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
    standard error (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall treatment effect and
    standard error (random effects model).}
  \item{seTE.predict}{Standard error utilised for prediction interval.}
  \item{lower.predict, upper.predict}{Lower and upper limits of prediction interval.}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance.}
  \item{method}{Pooling method: \code{"Inverse"}.}
  \item{call}{Function call.}
  \item{n.iter}{Actual number of iterations to estimate number of
    missing studies.}
  \item{trimfill}{A logical vector indicating studies that have been
    added by trim-and-fill method.}
  \item{df.hakn}{Degrees of freedom for test of treatment effect for
    Hartung-Knapp method (only if \code{hakn=TRUE}).}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{label.e}{Label for experimental group.}
  \item{label.c}{Label for control group.}
  \item{label.left}{Graph label on left side of forest plot.}
  \item{label.right}{Graph label on right side of forest plot.}
  \item{k0}{Number of studies added by trim-and-fill.}
  \item{n.e}{Number of observations in experimental group (only for
    object \code{x} of class \code{metabin} or \code{metacont}).}
  \item{n.c}{Number of observations in control group (only for object
    \code{x} of class \code{metabin} or \code{metacont}).}
  \item{event.e}{Number of events in experimental group (only for
    object \code{x} of class \code{metabin}).}
  \item{event.c}{Number of events in control group (only for object
    \code{x} of class \code{metabin}).}
  \item{mean.e}{Estimated mean in experimental group (only for object
    \code{x} of class \code{metacont}).}
  \item{sd.e}{Standard deviation in experimental group (only for
    object \code{x} of class \code{metacont}).}
  \item{mean.c}{Estimated mean in control group (only for object
    \code{x} of class \code{metacont}).}
  \item{sd.c}{Standard deviation in control group (only for object
    \code{x} of class \code{metacont}).}
  \item{n}{Number of observations (only for object \code{x} of class
    \code{metaprop}).}
  \item{event}{Number of events (only for object \code{x} of class
    \code{metaprop}).}
  \item{cor}{Corelation (only for object \code{x} of class
    \code{metacor}).}
  \item{class.x}{Main class of object \code{x} (e.g. 'metabin' or
    'metacont').}
  \item{version}{Version of R package meta used to create object.}
}

\references{
  Duval S & Tweedie R (2000a), A nonparametric "Trim and Fill" method of
  accounting for publication bias in meta-analysis. \emph{Journal of the
    American Statistical Association}, \bold{95}, 89--98.

  Duval S & Tweedie R (2000b), Trim and Fill: A simple funnel-plot-based
  method of testing and adjusting for publication bias in
  meta-analysis. \emph{Biometrics}, \bold{56}, 455--463.

  Peters JL, Sutton AJ, Jones DR, Abrams KR, Rushton L (2007),
  Performance of the trim and fill method in the presence of
  publication bias and between-study heterogeneity. \emph{Statisics in
  Medicine}, \bold{10}, 4544--62.
 
  Schwarzer G, Carpenter J, Rücker G (2010), Empirical evaluation
  suggests Copas selection model preferable to trim-and-fill method
  for selection bias in meta-analysis. \emph{Journal of Clinical
  Epidemiology}, \bold{63}, 282--8.

}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metagen}}, \code{\link{metabias}}, \code{\link{funnel}}}

\examples{
data(Fleiss93)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, sm="OR")
tf1 <- trimfill(meta1)
summary(tf1)
funnel(tf1)
funnel(tf1, pch=ifelse(tf1$trimfill, 1, 16),
       level=0.9, comb.random=FALSE)
#
# Use log odds ratios on x-axis
#
funnel(tf1, backtransf=FALSE)
funnel(tf1, pch=ifelse(tf1$trimfill, 1, 16),
       level=0.9, comb.random=FALSE, backtransf=FALSE)

trimfill(meta1$TE, meta1$seTE, sm=meta1$sm)
}

\keyword{htest}
