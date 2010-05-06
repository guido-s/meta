\name{forest}

\alias{forest}

\title{Generic function to produce a forest plot}

\description{
  Draws a forest plot in the active graphics window (using grid
  graphics system).
}

\usage{
forest(x, ...)
}
\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{\dots}{Additional arguments.}
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{forest.meta}}, \code{\link{metabin}},
  \code{\link{metacont}}, \code{\link{metagen}}, \code{\link[grid]{grid.newpage}} }

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", meth="I",
                 studlab=paste(author, year))

grid.newpage()
forest(meta1)
}
\keyword{hplot}
