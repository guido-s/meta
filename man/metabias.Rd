\name{metabias}
\alias{metabias}
\title{Generic function to test for funnel plot asymmetry}
\description{
  Test for funnel plot asymmetry, based on rank correlation or linear
  regression method.
}
\usage{
metabias(x, ...)
}
\arguments{
  \item{x}{An object of class \code{meta} or estimated treatment
    effect in individual studies.}
  \item{\dots}{Additional arguments.}
}
\details{
  For more details, see commands \code{\link{metabias.meta}} and
  \code{\link{metabias.default}}.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabias.meta}}, \code{\link{funnel}}, \code{\link{funnel.meta}}, \code{\link{metabin}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=1:10,
                 sm="RR", method="I")

##
## Using function metabias.meta:
##
metabias(meta1)

##
## Same result using function metabias.default:
##
metabias(meta1$TE, meta1$seTE)

}
\keyword{htest}
