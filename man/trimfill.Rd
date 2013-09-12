\name{trimfill}

\alias{trimfill}

\title{Generic function for trim-and-fill method}

\description{
  Trim and fill method for estimating and adjusting for the number and
  outcomes of missing studies in a meta-analysis.
}

\usage{
trimfill(x, ...)
}

\arguments{
  \item{x}{An object of class \code{meta}, or estimated treatment
    effect in individual studies.}
  \item{\dots}{Additional arguments as in \code{par}.}
}

\details{
  The trim and fill method (Duval, Tweedie 2000a, 2000b) can be used
  for estimating and adjusting for the number and outcomes of missing
  studies in a meta-analysis. The method relies on scrutiny of one
  side of a funnel plot for asymmetry assumed due to publication
  bias. For more details, see help page of R function
  \code{\link{trimfill.meta}}.
}

\value{
  An object of class \code{c("metagen", "meta", "trimfill")}. The object is a
  list containing the following components:
  \item{studlab, sm, left, ma.fixed, type}{}
  \item{n.iter.max, level, level.comb,}{As defined above.}
  \item{comb.fixed, comb.random}{}
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
  \item{k0}{Number of studies added by trim and fill.}
}

\references{
  Duval S & Tweedie R (2000a), A nonparametric "Trim and Fill" method of
  accounting for publication bias in meta-analysis. \emph{Journal of the
    American Statistical Association}, \bold{95}, 89--98.

  Duval S & Tweedie R (2000b), Trim and Fill: A simple funnel-plot-based
  method of testing and adjusting for publication bias in
  meta-analysis. \emph{Biometrics}, \bold{56}, 455--463.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metagen}}, \code{\link{metabias}}, \code{\link{trimfill.meta}}, \code{\link{funnel}}}

\examples{
data(Fleiss93)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, sm="OR")
tf1 <- trimfill(meta1)
summary(tf1)
funnel(tf1, pch=ifelse(tf1$trimfill, 1, 16),
       level=0.95, comb.fixed=TRUE)

trimfill(meta1$TE, meta1$seTE, sm=meta1$sm)
}

\keyword{htest}
