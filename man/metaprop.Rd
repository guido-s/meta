\name{metaprop}
\alias{metaprop}

\title{Meta-analysis of single proportions}

\description{
  Calculation of an overall proportion from studies reporting a single
  proportion.
}

\usage{
metaprop(event, n, studlab,
         data = NULL, subset = NULL,
         sm="PFT", freeman.tukey,
         incr=0.5, allincr=FALSE, addincr=FALSE,
         level = 0.95, level.comb = level,
         comb.fixed=TRUE, comb.random=TRUE,
         title="", complab="", outclab="",
         byvar, bylab, print.byvar=TRUE,
         warn=TRUE)
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
  \item{freeman.tukey}{Use of this argument is deprecated, please
    instead use argument \code{sm}, see Details.}
  \item{incr}{A numeric which is added to each cell frequency for
    studies with a zero cell count.}
  \item{allincr}{A logical indicating if \code{incr} is added to each
    cell frequency of all studies if at least one study has a zero cell
    count. If false, \code{incr} is added only to each cell frequency of
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
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{event.e}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{warn}{A logical indicating whether the addition of \code{incr}
    to studies with zero cell frequencies should result in a warning.}
}

\details{
  Fixed effect and random effects meta-analysis of single proportions to
  calculate an overall proportion.
  
  The following transformations of proportions are implemented to
  calculate an overall proportion:

  \itemize{
  \item \code{sm="PFT"}:  Freeman-Tukey Double arcsine transformation
  
  \item \code{sm="PAS"}:  Arcsine transformation
  
  \item \code{sm="PRAW"}:  Raw, i.e. untransformed, proportions
  
  \item \code{sm="PLN"}:  Log transformation
  
  \item \code{sm="PLOGIT"}:  Logit transformation}
  
  In older versions of the R package meta (< 1.5.0), only the
  Freeman-Tukey Double arcsine transformation and the arcsine
  transformation were implemenated and an argument \code{freeman.tukey}
  could be used to distinguish between these two methods. The use of the
  argument \code{freeman.tukey} is deprecated and the argument will
  probably removed from the function \code{metaprop} in the future. At
  the moment, a correponding warning message is printed if the argument
  \code{freeman.tukey} is used.
  
  If the summary measure is equal to "PRAW", "PLN", or "PLOGIT", a
  continuity correction is applied if any studies has a zero cell
  count. By default, 0.5 is added to all cell frequencies of studies
  with a zero cell count (argument \code{incr}).
  
  Internally, both fixed effect and random effects models are calculated
  regardless of values choosen for arguments \code{comb.fixed} and
  \code{comb.random}. Accordingly, the estimate for the random effects
  model can be extracted from component \code{TE.random} of an object
  of class \code{"meta"} even if \code{comb.random=FALSE}. However, all
  functions in R package \code{meta} will adequately consider the values
  for \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.meta}} will not print results for the random effects
  model if \code{comb.random=FALSE}.
}

\value{
  An object of class \code{c("metaprop", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{event, n, studlab,}{}
  \item{sm, incr, allincr, addincr, }{}
  \item{level, level.comb,}{As defined above.}
  \item{comb.fixed, comb.random,}{}
  \item{byvar, bylab, print.byvar, warn}{}
  \item{TE, seTE}{Estimated (un)transformed proportion and its standard
    error for individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall (un)transformed
    proportion and standard error (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall (un)transformed
    proportion and standard error (random effects model).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance (moment estimator of
    DerSimonian-Laird).}
  \item{sm}{A character string: \code{"proportion"}}
  \item{method}{A character string indicating method used
    for pooling: \code{"Inverse"}}
  \item{call}{Function call.}
  \item{version}{Version of R package meta used to create object.}
}

\references{
  Edward JM et al. (2006),
  Adherence to antiretroviral therapy in sub-saharan Africa and North
  America - a meta-analysis. \emph{Journal of the American Medical
  Association}, \bold{296}, 679--690.

  Freeman MF & Tukey JW (1950),
  Transformations related to the angular and the square
  root. \emph{Annals of Mathematical Statistics},
  \bold{21}, 607--611.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
metaprop(4:1, c(10, 20, 30, 40))
metaprop(4:1, c(10, 20, 30, 40), sm="PAS")
metaprop(4:1, c(10, 20, 30, 40), sm="PRAW")
metaprop(4:1, c(10, 20, 30, 40), sm="PLN")
metaprop(4:1, c(10, 20, 30, 40), sm="PLOGIT")
}

\keyword{htest}
