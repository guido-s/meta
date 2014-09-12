\name{read.mtv}
\alias{read.mtv}
\title{Import RevMan 4 data files (.mtv)}
\description{
  Reads a file created with RevMan 4 and creates a data frame from it.
}
\usage{
read.mtv(file)
}
\arguments{
  \item{file}{The name of a file to read data values from.}
}
\details{
  Reads a file created with RevMan 4 (Menu: "File" - "Export" -
  "Analysis data file...") and creates a data frame from it.
}
\value{
  A data frame containing the following components:
  \item{comp.no}{Comparison number.}
  \item{outcome.no}{Outcome number.}
  \item{group.no}{Group number.}
  \item{studlab}{Study label.}
  \item{year}{Year of publication.}
  \item{event.e}{Number of events in experimental group.}
  \item{n.e}{Number of observations in experimental group.}
  \item{event.c}{Number of events in control group.}
  \item{n.c}{Number of observations in control group.}
  \item{mean.e}{Estimated mean in experimental group.}
  \item{sd.e}{Standard deviation in experimental group.}
  \item{mean.c}{Estimated mean in control group.}
  \item{sd.c}{Standard deviation in control group.}
  \item{O.E}{Observed minus expected (IPD analysis).}
  \item{V}{Variance of \code{O.E} (IPD analysis).}
  \item{order}{Ordering of studies.}
  \item{conceal}{Concealment of treatment allocation.}
  \item{grplab}{Group label.}
  \item{type}{Type of outcome. D = dichotomous, C = continuous, P = IPD.}
  \item{outclab}{Outcome label.}
  \item{graph.exp}{Graph label for experimental group.}
  \item{graph.cont}{Graph label for control group.}
  \item{label.exp}{Label for experimental group.}
  \item{label.cont}{Label for control group.}
  \item{complab}{Comparison label.}
}
\references{
  \emph{Review Manager (RevMan)} [Computer program]. Version 4.2 for
  Windows. Copenhagen: The Nordic Cochrane Centre, The Cochrane
  Collaboration, 2003.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
# Locate MTV-data file "FLEISS93.MTV" in sub-directory of package "meta"
#
filename <- system.file("extdata/Fleiss93.MTV", package = "meta")
#
fleiss93.cc <- read.mtv(filename)

# Same result as R Command example(Fleiss93):
#
metabin(event.e, n.e, event.c, n.c,
        data=fleiss93.cc, subset=type=="D",
        studlab=paste(studlab, year))

# Same result: example(Fleiss93cont)
#
metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
        data=fleiss93.cc, subset=type=="C",
        studlab=paste(studlab, year))
}
\keyword{datagen}
