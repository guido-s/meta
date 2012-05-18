\name{trimfill.meta}

\alias{trimfill.meta}
\alias{trimfill.default}

\title{Trim and fill method for meta-analysis}

\description{
  Trim and fill method for estimating and adjusting for the number and
  outcomes of missing studies in a meta-analysis.
}

\usage{
\method{trimfill}{default}(x, seTE, left=NULL, ma.fixed=TRUE, type="L", n.iter.max=50,
         sm=NULL, studlab=NULL, level=0.95, level.comb=0.95,
         comb.fixed=TRUE, comb.random=TRUE,
         hakn=FALSE, method.tau="DL", silent=TRUE, ...)

\method{trimfill}{meta}(x, left=NULL, ma.fixed=TRUE, type="L", n.iter.max=50,
         sm=NULL, studlab=NULL, level=x$level, level.comb=x$level.comb,
         comb.fixed=x$comb.fixed, comb.random=x$comb.random,
         hakn=x$hakn, method.tau=x$method.tau, silent=TRUE, ...)
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
    measure, e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"AS"},
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
    to estimate the between-study variance tau-squared. Either
    \code{"DL"}, \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
    \code{"HE"}, or \code{"EB"}, can be abbreviated.}
  \item{silent}{A logical indicating whether basic information on
    iterations shown.}
  \item{...}{other arguments}
}

\details{
  The trim and fill method can be used for estimating and adjusting for
  the number and outcomes of missing studies in a meta-analysis. The
  method relies on scrutiny of one side of a funnel plot for asymmetry
  assumed due to publication bias.

  Internally, both fixed effect and random effects models are calculated
  regardless of values choosen for arguments \code{comb.fixed} and
  \code{comb.random}. Accordingly, the estimate for the random effects
  model can be extracted from component \code{TE.random} of an object
  of class \code{"meta"} even if \code{comb.random=FALSE}. However, all
  functions in R package \code{meta} will adequately consider the values
  for \code{comb.fixed} and \code{comb.random}. E.g. function
  \code{\link{print.meta}} will not print results for the random effects
  model if \code{comb.random=FALSE}.

  The function \code{metagen} is called internally.
}

\value{
  An object of class \code{c("metagen", "meta", "trimfill")}. The object is a
  list containing the following components:
  \item{studlab, sm, left, ma.fixed, type,}{}
  \item{n.iter.max, level, level.comb,}{As defined above.}
  \item{comb.fixed, comb.random, hakn, method.tau}{}
  \item{TE, seTE}{Estimated treatment effect and standard error of individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
    standard error (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall treatment effect and
    standard error (random effects model).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance.}
  \item{method}{Pooling method: \code{"Inverse"}.}
  \item{call}{Function call.}
  \item{n.iter}{Actual number of iterations to estimate number of
    missing studies.}
  \item{trimfill}{A logical vector indicating studies that have been
    added by trim and fill method.}
  \item{df.hakn}{Degrees of freedom for test of treatment effect for
    Hartung-Knapp method (only if \code{hakn=TRUE}).}
  \item{k0}{Number of studies added by trim and fill.}
  \item{version}{Version of R package meta used to create object.}
}

\references{
  Duval S & Tweedie R (2000), A nonparametric "Trim and Fill" method of
  accounting for publication bias in meta-analysis. \emph{Journal of the
    American Statistical Association}, \bold{95}, 89--98.

  Duval S & Tweedie R (2000), Trim and Fill: A simple funnel-plot-based
  method of testing and adjusting for publication bias in
  meta-analysis. \emph{Biometrics}, \bold{56}, 455--463.
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

trimfill(meta1$TE, meta1$seTE, sm=meta1$sm)
}

\keyword{htest}
