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
         data=NULL, subset=NULL, sm="MD",
         level = 0.95, level.comb = level,
         comb.fixed=TRUE, comb.random=TRUE,
         title="", complab="", outclab="",
         label.e="Experimental", label.c="Control",
         byvar, bylab, print.byvar=TRUE)
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
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{label.e}{Label for experimental group.}
  \item{label.c}{Label for control group.}
  \item{sm}{A character string indicating which summary measure
    (\code{"MD"} or \code{"SMD"}) is to be used for pooling of
    studies.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{n.e}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
}

\details{
  Calculation of fixed and random effects estimates for meta-analyses
  with continuous outcome data; inverse variance weighting is used for
  pooling. The DerSimonian-Laird estimate is used in the
  random effects model. The mean difference is used as measure of
  treatment effect if \code{sm="MD"} -- which correspond to
  \code{sm="WMD"} in older versions (<0.9) of the meta package. For the
  summary measure \code{"SMD"}, Hedges' adjusted g is utilised for
  pooling.

  Internally, both fixed effect and random effects models are calculated
  regardless of values choosen for arguments \code{comb.fixed} and
  \code{comb.random}. Accordingly, the estimate for the random effects
  model can be extracted from component \code{TE.random} of an object
  of class \code{"meta"} even if \code{comb.random=FALSE}. However, all
  functions in R package \code{meta} will adequately consider the values
  for \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.meta}} will not print results for the random effects
  model if \code{comb.random=FALSE}.

  The function \code{metagen} is called internally to calculate
  individual and overall treatment estimates and standard errors.  
}

\value{
  An object of class \code{c("metacont", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{n.e, mean.e, sd.e,}{}
  \item{n.c, mean.c, sd.c,}{}
  \item{studlab, sm, level, level.comb,}{}
  \item{comb.fixed, comb.random,}{}
  \item{byvar, bylab, print.byvar}{As defined above.}
  \item{TE, seTE}{Estimated treatment effect and standard error of individual studies.}
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
  \item{version}{Version of R package meta used to create object.}
}

\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metagen}}}

\examples{
data(Fleiss93cont)
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont, sm="SMD")
meta1

meta2 <- metacont(Fleiss93cont$n.e, Fleiss93cont$mean.e,
                  Fleiss93cont$sd.e,
                  Fleiss93cont$n.c, Fleiss93cont$mean.c,
                  Fleiss93cont$sd.c,
                  sm="SMD")
meta2
}

\keyword{htest}
