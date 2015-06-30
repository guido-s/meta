\name{metacr}

\alias{metacr}

\title{Meta-analysis of outcome data from Cochrane review}

\description{
  Wrapper function to perform meta-analysis for a single outcome of a
  Cochrane Intervention review.
}

\usage{
metacr(x, comp.no=1, outcome.no=1,
       method, sm,
       level=.settings$level, level.comb=.settings$level.comb,
       comb.fixed, comb.random,
       hakn=FALSE,
       method.tau="DL",
       tau.common=FALSE,
       prediction=.settings$prediction, level.predict=.settings$level.predict,
       swap.events, logscale,
       backtransf=.settings$backtransf,
       title, complab, outclab,
       keepdata=.settings$keepdata, warn=FALSE)
}

\arguments{
  \item{x}{An object of class \code{rm5} created by R function
  \code{read.rm5}.}
  \item{comp.no}{Comparison number.}
  \item{outcome.no}{Outcome number.}
  \item{method}{A character string indicating which method is to be used
    for pooling of studies. One of \code{"Inverse"}, \code{"MH"}, or
    \code{"Peto"}, can be abbreviated.}
  \item{sm}{A character string indicating which summary measure
    (\code{"RR"}, \code{"OR"}, \code{"RD"}, \code{"ASD"}, \code{"HR"},
    \code{"MD"}, or \code{"SMD"}) is to be used for pooling of
    studies.}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether a fixed effect
    meta-analysis should be conducted.}
  \item{comb.random}{A logical indicating whether a random effects
    meta-analysis should be conducted.}
  \item{hakn}{A logical indicating whether the method by Hartung and
          Knapp should be used to adjust test statistics and
          confidence intervals.}
  \item{method.tau}{A character string indicating which method is used
    to estimate the between-study variance \eqn{\tau^2}. Either
    \code{"DL"}, \code{"PM"}, \code{"REML"}, \code{"ML"}, \code{"HS"},
    \code{"SJ"}, \code{"HE"}, or \code{"EB"}, can be abbreviated.}
  \item{tau.common}{A logical indicating whether tau-squared should be
    the same across subgroups.}
  \item{prediction}{A logical indicating whether a prediction interval
    should be printed.}
  \item{level.predict}{The level used to calculate prediction interval
    for a new study.}
  \item{swap.events}{A logical indicating whether events and
    non-events should be interchanged.}
  \item{logscale}{A logical indicating whether effect estimates are
    entered on log-scale.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and plots. If \code{backtransf=TRUE}
    (default), results for \code{sm="OR"} are printed as odds ratios
    rather than log odds ratios and results for \code{sm="ZCOR"} are
    printed as correlations rather than Fisher's z transformed
    correlations, for example.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{keepdata}{A logical indicating whether original data (set)
    should be kept in meta object.}
  \item{warn}{A logical indicating whether warnings should be printed
    (e.g., if \code{incr} is added to studies with zero cell
    frequencies).}
}

\details{
  Cochrane Intervention reviews are based on the comparison of two
  interventions. Each Cochrane Intervention review can have a variable
  number of comparisons. For each comparison, a variable number of
  outcomes can be define. For each outcome, a seperate meta-analysis is
  conducted. Review Manager 5 (RevMan 5) is the current software used for
  preparing and maintaining Cochrane Reviews
  (\url{http://tech.cochrane.org/revman/}).
  
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
  \emph{Review Manager (RevMan)} [Computer program]. Version
  5.3. Copenhagen: The Nordic Cochrane Centre, The Cochrane
  Collaboration, 2014.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{read.rm5}}}

\examples{
# Locate export data file "Fleiss93_CR.csv"
# in sub-directory of package "meta"
#
filename <- system.file("data/Fleiss93_CR.csv.gz", package = "meta")
#
Fleiss93_CR <- read.rm5(filename)

# Same result as R command example(Fleiss93):
#
metacr(Fleiss93_CR)

# Same result as R command example(Fleiss93cont):
#
metacr(Fleiss93_CR, 1, 2)

forest(metacr(Fleiss93_CR, 1, 2))

# Change summary measure to RR
#
m1 <- metacr(Fleiss93_CR)
update(m1, sm="RR")

}

\keyword{htest}
