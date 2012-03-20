\name{ci}
\alias{ci}
\title{Calculation of confidence intervals (based on normal
approximation or t-distribution)}
\description{
  Calculation of confidence intervals; based on normal approximation
  or t-distribution.
}
\usage{
ci(TE, seTE, level=0.95, df=NULL)
}
\arguments{
  \item{TE}{Estimated treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{level}{The confidence level required.}
  \item{df}{Degrees of freedom (for confidence intervals based on t-distribution).}
}

\value{
  List with components
  \item{TE}{Estimated treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{lower}{Lower confidence limits.}
  \item{upper}{Upper confidence limits.}  
  \item{z}{Test statistic (either z-score or t-score).}
  \item{p}{P-value of test with null hypothesis \code{TE=0}.}
  \item{level}{The confidence level required.}
  \item{df}{Degrees of freedom (t-distribution).}
}

\note{
  This function is primarily called from other functions of the library
  \code{meta}, e.g. \code{forest.meta}, \code{summary.meta}.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\examples{
data.frame(ci(170, 10))
data.frame(ci(170, 10, 0.99))
data.frame(ci(1.959964, 1))
data.frame(ci(2.2621571628, 1, df=9))
}
\keyword{htest}
