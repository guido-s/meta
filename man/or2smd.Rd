% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/or2smd.R
\name{or2smd}
\alias{or2smd}
\title{Conversion from log odds ratio to standardised mean difference}
\usage{
or2smd(
  lnOR,
  selnOR,
  studlab,
  data = NULL,
  subset = NULL,
  exclude = NULL,
  method = "HH",
  ...
)
}
\arguments{
\item{lnOR}{Log odds ratio(s) or meta-analysis object.}

\item{selnOR}{Standard error(s) of log odds ratio(s) (ignored if
argument \code{lnOR} is a meta-analysis object).}

\item{studlab}{An optional vector with study labels (ignored if
argument \code{lnOR} is a meta-analysis object).}

\item{data}{An optional data frame containing the study information
(ignored if argument \code{lnOR} is a meta-analysis object).}

\item{subset}{An optional vector specifying a subset of studies to
be used (ignored if argument \code{lnOR} is a meta-analysis
object).}

\item{exclude}{An optional vector specifying studies to exclude
from meta-analysis, however, to include in printouts and forest
plots (ignored if argument \code{lnOR} is a meta-analysis
object).}

\item{method}{A character string indicating which method is used to
convert log odds ratios to standardised mean differences. Either
\code{"HH"} or \code{"CS"}, can be abbreviated.}

\item{\dots}{Additional arguments passed on to
\code{\link{metagen}} (ignored if argument \code{lnOR} is a
meta-analysis object).}
}
\value{
An object of class \code{c("metagen", "meta")} with corresponding
generic functions (see \code{\link{meta-object}}).
}
\description{
Conversion from log odds ratio to standardised mean difference
using method by Hasselblad & Hedges (1995) or Cox (1970).
}
\details{
This function implements the following methods for the conversion
from log odds ratios to standardised mean difference:
\itemize{
\item Hasselblad & Hedges (1995) assuming logistic distributions
  (\code{method == "HH"})
\item Cox (1970) and Cox & Snell (1989) assuming normal
  distributions (\code{method == "CS"})
}
Internally, \code{\link{metagen}} is used to conduct a
meta-analysis with the standardised mean difference as summary
measure.

Argument \code{lnOR} can be either a vector of log odds ratios or a
meta-analysis object created with \code{\link{metabin}} or
\code{\link{metagen}} and the odds ratio as summary measure.

Argument \code{selnOR} is mandatory if argument \code{lnOR} is a
vector and ignored otherwise. Additional arguments in \code{\dots}
are only passed on to \code{\link{metagen}} if argument \code{lnOR}
is a vector.
}
\examples{
# Example from Borenstein et al. (2009), Chapter 7
#
mb <- or2smd(0.9069, sqrt(0.0676))
# TE = standardised mean difference (SMD); seTE = standard error of SMD
data.frame(SMD = round(mb$TE, 4), varSMD = round(mb$seTE^2, 4))

# Use dataset from Fleiss (1993)
#
data(Fleiss1993bin)
m1 <- metabin(d.asp, n.asp, d.plac, n.plac,
  data = Fleiss1993bin, studlab = paste(study, year),
  sm = "OR", random = FALSE)
or2smd(m1)

}
\references{
Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009):
\emph{Introduction to Meta-Analysis}.
Chichester: Wiley

Cox DR (1970):
\emph{Analysis of Binary Data}.
London: Chapman and Hall / CRC

Cox DR, Snell EJ (1989):
\emph{Analysis of Binary Data} (2nd edition).
London: Chapman and Hall / CRC

Hasselblad V, Hedges LV (1995):
Meta-analysis of screening and diagnostic tests.
\emph{Psychological Bulletin},
\bold{117}, 167--78
}
\seealso{
\code{\link{smd2or}}, \code{\link{metabin}},
  \code{\link{metagen}}, \code{\link{metacont}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
