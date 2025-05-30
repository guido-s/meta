% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bracket.longarm.R
\name{[.longarm}
\alias{[.longarm}
\alias{[,longarm}
\title{Extract parts of longarm object}
\usage{
\method{[}{longarm}(x, ...)
}
\arguments{
\item{x}{An object of class \code{\link{longarm}}.}

\item{\dots}{Additional arguments (passed on to [.data.frame).}
}
\description{
Auxiliary function to extract parts of \code{\link{longarm}} object.
}
\examples{
# Transform data from wide arm-based format to contrast-based format
pw1 <- pairwise(list(Treatment1, Treatment2, Treatment3),
  n = list(n1, n2, n3),
  mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3),
  data = dat.franchini2012, studlab = Study)

# Transform data from contrast-based to long arm-based format
# and only keep the main variables
la1 <- longarm(pw1, append = FALSE)
head(la1)

la1[la1$studlab == "Lieberman 1998", ]

}
\seealso{
\code{\link{longarm}}, \code{\link{subset.longarm}},
  \code{\link[metadat]{dat.franchini2012}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
