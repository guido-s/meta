\name{metacor}
\alias{metacor}

\title{Meta-analysis of correlations}

\description{
  Calculation of fixed and random effects estimates for meta-analyses
  with correlations; inverse variance weighting is used for pooling.
}

\usage{
metacor(cor, n, studlab,
        data=NULL, subset=NULL,
        sm="ZCOR",
        level=0.95, level.comb=level,
        comb.fixed=TRUE, comb.random=TRUE,
        title="", complab="", outclab="",
        byvar, bylab, print.byvar=TRUE
        )
}

\arguments{
  \item{cor}{Correlation.}
  \item{n}{Number of observations.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., cor and n.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating which summary measure
    (\code{"ZCOR"} or \code{"COR"}) is to be used for pooling of
    studies.}
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
  Fixed effect and random effects meta-analysis of correlations
  based either on Fisher's z transformation of proportions (\code{sm="ZCOR"})
  or direct combination of correlations (\code{sm="COR"}) (see Cooper et
  al., p264-5 and p273-4).

  Only few statisticians would advocate the use of untransformed
  correlations unless sample sizes are very large (see Cooper et al.,
  p265). The artificial example given below shows that the smallest
  study gets the largest weight if correlations are combined directly
  because the correlation is closest to 1.

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
  An object of class \code{c("metacor", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{cor, n, studlab,}{}
  \item{sm, level, level.comb,}{}
  \item{comb.fixed, comb.random,}{}
  \item{byvar, bylab, print.byvar}{As defined above.}
  \item{TE, seTE}{Either Fisher's z transformation of correlations
    (\code{sm="ZCOR"}) or correlations (\code{sm="COR"}) for individual
    studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall effect (Fisher's z
    transformation of correlation or correlation) and its standard error
    (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall effect (Fisher's z
    transformation of correlation or correlation) and its standard error
    (random effects model).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance (moment estimator of
    DerSimonian-Laird).}
  \item{method}{A character string indicating method used
    for pooling: \code{"Inverse"}}
  \item{call}{Function call.}
  \item{version}{Version of R package meta used to create object.}
}

\references{
  Cooper H, Hedges LV, Valentine JC (2009),
  \emph{The Handbook of Research Synthesis and Meta-Analysis}, 2nd Edition.
  New York: Russell Sage Foundation.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metacont}}, \code{\link{metagen}}, \code{\link{print.meta}}}

\examples{
metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))
metacor(c(0.85, 0.7, 0.95), c(20, 40, 10), sm="cor")
}

\keyword{htest}
