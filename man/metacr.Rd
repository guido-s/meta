\name{metacr}

\alias{metacr}

\title{Meta-analysis of outcome data from Cochrane review}

\description{
  Wrapper function to perform meta-analysis for a single outcome of a
  Cochrane Intervention review.
}

\usage{
metacr(x, comp.no=1, outcome.no=1, smother="", logscale=TRUE)
}

\arguments{
  \item{x}{An object of class \code{rm5} created by R function
  \code{read.rm5}.}
  \item{comp.no}{Comparison number.}
  \item{outcome.no}{Outcome number.}
  \item{smother}{A character specifying the summary measure to use
    instead of the value "OTHER"}
  \item{logscale}{A logical indicating whether treatment effects have
    been entered on the log scale in RevMan5. Only relevant for an
    outcome with summary measure equal to "OTHER".}
}

\details{
  Cochrane Intervention reviews are based on the comparison of two
  interventions. Each Cochrane Intervention review can have a variable
  number of comparisons. For each comparison, a variable number of
  outcomes can be define. For each outcome, a seperate meta-analysis is
  conducted. Review Manager 5 (RevMan 5) is the current software used for
  preparing and maintaining Cochrane Reviews
  (\url{http://www.cc-ims.net/revman/}).
  
  This wrapper function can be used to perform meta-analysis for a
  single outcome of a Cochrane Intervention review. Internally, R
  functions \code{metabin}, \code{metacont}, and \code{metagen} are
  called - depending on the definition of the outcome in RevMan 5.
}

\value{
  An object of class \code{"meta"} and \code{"metabin"},
  \code{"metacont"}, or \code{"metagen"} depending on outcome type
  utilised in Cochrane Intervention review for selected outcome.
}

\references{
  \emph{Review Manager (RevMan)} [Computer program]. Version 5.0, The Nordic
  Cochrane Centre, The Cochrane Collaboration, 2008.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{read.rm5}}}

\examples{
## Locate export data file "Fleiss93_CR.csv" or "Fleiss93_CR_Windows.csv"
## in sub-directory of package "meta"
##
filename <- paste(searchpaths()[seq(along=search())[search()==
                  "package:meta"]], "/data/Fleiss93_CR",
                  if (Sys.info()[["sysname"]]=="Windows") "_Windows" else "",
                  ".csv", sep="")
Fleiss93_CR <- read.rm5(filename)

## Same result: example(Fleiss93)
##
metacr(Fleiss93_CR)

## Same result: example(Fleiss93cont)
##
metacr(Fleiss93_CR, 1, 2)
}

\keyword{htest}
