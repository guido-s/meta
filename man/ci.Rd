\name{ci}
\alias{ci}
\title{Calculation of confidence intervals (normal approximation)}
\description{
  Calculation of confidence intervals; based on normal approximation.}
\usage{
ci(TE, seTE, level=0.95)
}
\arguments{
  \item{TE}{Estimated treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{level}{The confidence level required.}
}

\value{
  List with components
  \item{TE}{Estimated treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{lower}{Lower confidence limits.}
  \item{upper}{Upper confidence limits.}  
  \item{zscore}{Test statistic.}
  \item{p}{P-value of test with null hypothesis \code{TE=0}.}
  \item{level}{The confidence level required.}
}

\note{
  This function is primarily called from other functions of the library
  \code{meta}, e.g. \code{plot.meta}, \code{summary.meta}.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\examples{
as.data.frame(ci(170, 10))
as.data.frame(ci(170, 10, 0.99))
}
\keyword{htest}

