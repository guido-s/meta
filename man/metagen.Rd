\name{metagen}

\alias{metagen}

\title{Generic inverse variance meta-analysis}

\description{
  Fixed and random effects meta-analysis based on estimates (e.g. log
  hazard ratios) and their standard errors; inverse variance weighting
  is used for pooling.
}

\usage{
metagen(TE, seTE, studlab, data=NULL, subset=NULL, sm="",
        level = 0.95, level.comb = level,
        comb.fixed=TRUE, comb.random=TRUE,
        title="", complab="", outclab="",
        label.e="Experimental", label.c="Control",
        byvar, bylab, print.byvar=TRUE)
}

\arguments{
  \item{TE}{Estimate of treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"AS"},
    \code{"MD"}, \code{"SMD"}.}
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
  \item{label.e}{Label for experimental group.}
  \item{label.c}{Label for control group.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{TE}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
}

\details{
  Generic method for meta-analysis, only treatment estimates and their
  standard error are needed. The method is useful, e.g., for pooling of
  survival data (using log hazard ratio and standard errors as
  input). The inverse variance method is used for pooling. Random
  effects estimate is based on the DerSimonian-Laird method.

  Internally, both fixed effect and random effects models are
  calculated regardless of values choosen for arguments
  \code{comb.fixed} and \code{comb.random}. Accordingly, the estimate
  for the random effects model can be extracted from component
  \code{TE.random} of an object of class \code{"meta"} even if
  \code{comb.random=FALSE}. However, all functions in R package
  \code{meta} will adequately consider the values for \code{comb.fixed}
  and \code{comb.random}. E.g. function \code{\link{print.meta}} will
  not print results for the random effects model if
  \code{comb.random=FALSE}.
}

\value{
  An object of class \code{c("metagen", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{TE, seTE, studlab,}{}
  \item{sm, level, level.comb,}{}
  \item{comb.fixed, comb.random,}{}
  \item{byvar, bylab, print.byvar}{As defined above.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
    standard error (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall treatment effect and
    standard error (random effects model).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic.}
  \item{tau}{Square-root of between-study variance (moment estimator of
    DerSimonian-Laird).}
  \item{method}{Pooling method: \code{"Inverse"}.}
  \item{call}{Function call.}
}

\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{print.meta}}}

\examples{
data(Fleiss93)
meta1 <- metabin(event.e, n.e, event.c, n.c, data=Fleiss93, sm="RR", meth="I")
meta1

##
## Identical results by using the following commands:
##
meta1
metagen(meta1$TE, meta1$seTE, sm="RR")


##
## Meta-analysis of survival data:
##
logHR <- log(c(0.95, 1.5))
selogHR <- c(0.25, 0.35)

metagen(logHR, selogHR, sm="HR")
}

\keyword{htest}
