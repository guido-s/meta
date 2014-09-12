\name{settings.meta}
\alias{settings.meta}
\title{Print and change default settings for meta-analyses in R package meta.}
\description{
  Print and change default settings for meta-analyses in R package meta.
}
\usage{
settings.meta(...)
}
\arguments{
  \item{...}{Arguments to change default settings.}
}

\details{

  This function can be used to define defaults for several arguments
  of the following R functions: \code{\link{metabin}},
  \code{\link{metacont}}, \code{\link{metacor}}, \code{\link{metacr}},
  \code{\link{metagen}}, \code{\link{metainc}},
  \code{\link{metaprop}}.
  
  A list of all arguments with current settings is printed using the
  command \code{settings.meta(print=TRUE)}.

  In order to reset all settings to the default settings of R package
  meta the command \code{settings.meta(reset=TRUE)} can be used.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\examples{
#
# Get list of settings
#
settings.meta(print=TRUE)

#
# Meta-analyses using default settings of R pacakge meta
#
metabin(10, 20, 15, 20, warn=FALSE)
metaprop(4, 20)
metabin(10, 20, 15, 20, sm="RD", warn=FALSE)
metaprop(4, 20, sm="PLN")
#
# Change default summary measure for R functions metabin and metaprop
#
oldset <- settings.meta(smbin="RD", smprop="PLN")
#
metabin(10, 20, 15, 20, warn=FALSE)
metaprop(4, 20)
#
# Reset settings
#
settings.meta(oldset)

#
# Change level for confidence intervals
#
metagen(1:3, (2:4)/10, sm="MD")
settings.meta(level=0.99, level.comb=0.999)
metagen(1:3, (2:4)/10, sm="MD")

#
# Always print a prediction interval
#
settings.meta(prediction=TRUE)
metagen(1:3, (2:4)/10, sm="MD")

#
# Try to set unknown argument
#
settings.meta(unknownarg=TRUE)

#
# Reset to default settings of R package meta
#
settings.meta(reset=TRUE)
metabin(10, 20, 15, 20, warn=FALSE)
metaprop(4, 20)
metagen(1:3, (2:4)/10, sm="MD")

#
# Never back transform results (e.g. print log odds ratios instead of
# odds ratios, print transformed correlations/proportions instead of
# correlations/proportions)
#
oldset <- settings.meta(backtransf=FALSE)
metabin(10, 20, 15, 20, warn=FALSE)
metaprop(4, 20)
metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))
settings.meta(oldset)
}
