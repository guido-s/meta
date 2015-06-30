\name{print.rm5}

\alias{print.rm5}
\alias{print.summary.rm5}
\alias{summary.rm5}
\alias{metabias.rm5}

\title{Print and summary methods for objects of class rm5}

\description{
  Print and summary methods for objects of class \code{rm5}.
}

\usage{
\method{print}{rm5}(x, ...)

\method{print}{summary.rm5}(x, ...)

\method{summary}{rm5}(object, comp.no, outcome.no, ...)

\method{metabias}{rm5}(x, comp.no, outcome.no,
                       method.bias="linreg",
                       method.bias.binary=method.bias,
                       method.bias.or="score",
                       k.min=10, ...)
}

\arguments{
  \item{x}{An object of class \code{rm5}.}
  \item{object}{An object of class \code{rm5}.}
  \item{comp.no}{Comparison number.}
  \item{outcome.no}{Outcome number.}
  \item{method.bias}{A character string indicating which test for
    small-study effects is to be used for all outcomes. Either
    \code{"rank"}, \code{"linreg"}, or \code{"mm"}, can be
    abbreviated. See function \code{\link{metabias}}}
  \item{method.bias.binary}{A character string indicating which test
    is to be used for binary outcomes. Either \code{"rank"},
    \code{"linreg"}, \code{"mm"}, \code{"count"}, \code{"score"}, or
    \code{"peters"}, can be abbreviated. See function
    \code{\link{metabias}}}
  \item{method.bias.or}{A character string indicating which test is to
    be used for binary outcomes with odds ratio as summary
    measure. Either \code{"rank"}, \code{"linreg"}, \code{"mm"},
    \code{"count"}, \code{"score"}, or \code{"peters"}, can be
    abbreviated. See function \code{\link{metabias}}}
  \item{k.min}{Minimum number of studies to perform test for small-study effects.}
  \item{...}{Additional arguments (ignored at the moment)}
}

\details{
  Review Manager 5 (RevMan 5) is the current software used for
  preparing and maintaining Cochrane Reviews
  (\url{http://tech.cochrane.org/revman/}). In RevMan 5, subgroup analyses
  can be defined and data from a Cochrane review can be imported to R
  using the function \code{read.rm5}.

  The R function \code{summary.rm5} can be used to redo all
  meta-analyses of the imported Cochrane Review.

  The R function \code{metabias.rm5} can be used to conduct a test for
  funnel plot asymmetry for all meta-analyses of the imported Cochrane
  Review.

  The R function \code{metacr} is called internally.
}

\references{
  Higgins, J.P.T and S. Green (2011),
  \emph{Cochrane Handbook for Systematic Reviews of Interventions Version 5.1.0 [Updated March 2011]}.
  The Cochrane Library: http://www.cochrane-handbook.org
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabias.meta}}, \code{\link{summary.meta}}, \code{\link{read.rm5}}}

\examples{
# Locate export data file "Fleiss93_CR.csv"
# in sub-directory of package "meta"
#
filename <- system.file("data/Fleiss93_CR.csv.gz", package = "meta")
#
Fleiss93_CR <- read.rm5(filename)

#
# Print summary results for all meta-analysis:
#
summary(Fleiss93_CR)

#
# Print results for tests of small-study effects:
#
metabias(Fleiss93_CR, k.min=5)
}

\keyword{print}
