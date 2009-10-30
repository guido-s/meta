\name{labbe}
\alias{labbe}
\title{L\'Abbe plot}
\description{
  Generic function for drawing a L\'Abbe plot.
}
\usage{
labbe(x, y, ...)
}
\arguments{
  \item{x}{The x coordinates of points of the L\'Abbe
    plot. Alternatively, an object of class \code{metabin}.}
  \item{y}{The y coordinates of the L\'Abbe plot, optional if \code{x}
    is an appropriate structure.}
  \item{\dots}{Parameters used in other L\'Abbe plot functions.}
}
\details{
  Generic function for drawing a L\'Abbe plot.
}
\references{
  L'Abbe  KA, Detsky AS, O'Rourke K (1987),
  Meta-analysis in clinical research. \emph{Annals of Internal
  Medicine}, \bold{107}, 224--233.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{labbe.metabin}}, \code{\link{metabin}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95,
                 studlab=paste(author, year),
                 sm="RR")


##
## L'Abbe plot
##
labbe(meta1)
}
\keyword{hplot}
