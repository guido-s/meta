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
         freeman.tukey=TRUE,
         level = 0.95, level.comb = level,
         comb.fixed=TRUE, comb.random=TRUE,
         title="", complab="", outclab="",
         byvar, bylab, print.byvar=TRUE)
}

\arguments{
  \item{event}{Number of events.}
  \item{n}{Number of observations.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., event and n.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{freeman.tukey}{A logical indicating if the Freeman-Tukey Double
    arcsine transformation should be used; otherwise the arcsine
    transformation is used.}
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
}

\details{
  Fixed effect and random effects meta-analysis of single proportions
  using either the Freeman-Tukey Double arcsine transformation or the
  arcsine transformation of proportions to calculate an overall
  proportion.

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
  \item{freeman.tukey, level, level.comb,}{}
  \item{comb.fixed, comb.random,}{}
  \item{byvar, bylab, print.byvar}{As defined above.}
  \item{TE, seTE}{Arcsine transformation of proportion and its standard
  error for individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall arcsine transformed
  proportion and standard error (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall arcsine transformed
  proportion and standard error (random effects model).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance (moment estimator of
    DerSimonian-Laird).}
  \item{sm}{A character string: \code{"proportion"}}
  \item{method}{A character string indicating method used
    for pooling: \code{"Inverse"}}
  \item{call}{Function call.}
}

\references{
  Freeman MF & Tukey JW (1950),
  Transformations related to the angular and the square
  root. \emph{Annals of Mathematical Statistics},
  \bold{21}, 607--611.

  Miller JJ (1978),
  The inverse of the Freeman-Tukey double arcsine
  transformation. \emph{The American Statistician},
  \bold{32}, 138.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
metaprop(0, 10)
metaprop(0, 10, freeman.tukey=FALSE)
}

\keyword{htest}
