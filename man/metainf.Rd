\name{metainf}
\alias{metainf}
\title{Influence analysis in meta-analysis}
\description{
  Performs a influence analysis. Pooled estimates are calculated
  omitting one study at a time.
}
\usage{
metainf(x, pooled, sortvar)
}
\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{pooled}{A character string indicating whether a fixed effect or
  random effects model is used for pooling. Either missing (see
  Details), \code{"fixed"} or \code{"random"}, can be abbreviated.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
}
\details{
  Performs a influence analysis; pooled estimates are calculated
  omitting one study at a time. Studies are sorted according to \code{sortvar}.

  Information from object \code{x} is utilised if argument \code{pooled}
  is missing. A fixed effect model is assumed (\code{pooled="fixed"}) if
  argument \code{x$comb.fixed} is \code{TRUE}; a random effects model
  is assumed (\code{pooled="random"}) if
  argument \code{x$comb.random} is \code{TRUE} and \code{x$comb.fixed}
  is \code{FALSE}.
}
\value{
  An object of class \code{c("metainf", "meta")} with corresponding
  \code{print}, \code{plot} function. The object is a
  list containing the following components:
  \item{TE, seTE}{Estimated treatment effect and standard error of
    pooled estimate in influence analysis.}
  \item{lower, upper}{Lower and upper confidence interval limits.}
  \item{studlab}{Study label describing omission of studies.}
  \item{p.value}{P-value for test of overall effect.}
  \item{w}{Sum of weights from fixed effect or random effects model.}
  \item{I2}{Heterogeneity statistic I2.}
  \item{tau}{Square-root of between-study variance.}
  \item{df.hakn}{Degrees of freedom for test of treatment effect for
  Hartung-Knapp method (only if \code{hakn=TRUE}).}
  \item{sm}{Summary measure.}
  \item{method}{Method used for pooling.}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{pooled}{As defined above.}
  \item{comb.fixed}{A logical indicating whether analysis is based on
    fixed effect model.}
  \item{comb.random}{A logical indicating whether analysis is based on
    random effects model.}
  \item{TE.fixed, seTE.fixed}{Value is \code{NA}.}
  \item{TE.random, seTE.random}{Value is \code{NA}.}
  \item{Q}{Value is \code{NA}.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{hakn}{A logical indicating whether the method by Hartung and
          Knapp is used to adjust test statistics and
          confidence intervals.}
  \item{method.tau}{A character string indicating which method is used
    to estimate the between-study variance \eqn{\tau^2}.}
  \item{tau.preset}{Prespecified value for the square-root of the
    between-study variance \eqn{\tau^2}.}
  \item{TE.tau}{Overall treatment effect used to estimate the
  between-study variance \eqn{\tau^2}.}
  \item{n.harmonic.mean}{Harmonic mean of number of observations (for back transformation of Freeman-Tukey Double arcsine transformation).}
  \item{version}{Version of R package meta used to create object.}
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
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, studlab=study,
                 sm="RR", method="I")
meta1

metainf(meta1)
metainf(meta1, pooled="random")

forest(metainf(meta1))
forest(metainf(meta1), layout="revman5")
forest(metainf(meta1, pooled="random"))

metainf(meta1, sortvar=study)
metainf(meta1, sortvar=7:1)

meta2 <- update(meta1, title="Fleiss93 meta-analysis",
                backtransf=FALSE)
metainf(meta2)
}
\keyword{htest}
