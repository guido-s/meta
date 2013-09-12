\name{addvar}
\alias{addvar}
\alias{as.data.frame.meta}
\title{Additional functions for objects of class meta}
\description{
  The \code{as.data.frame} method returns a data frame containing
  information on individual studies, e.g., estimated treatment effect
  and its standard error. The function \code{addvar} can be used to add
  a single variable to an object of class \code{meta} which for example
  is useful to add columns to a forest plot.
}
\usage{
\method{as.data.frame}{meta}(x, row.names=NULL, optional=FALSE, ...)

addvar(x, y, varname, by.x="studlab", by.y=by.x)
}
\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{row.names}{\code{NULL} or a character vector giving the row
  names for the data frame.}
  \item{optional}{logical. If \code{TRUE}, setting row names and
    converting column names (to syntactic names) is optional.}
  \item{y}{A data frame with an additional covariate}
  \item{varname}{A character specifying name of additional variable}
  \item{by.x, by.y}{Specifications of the common columns (see \code{merge})}
  \item{\dots}{other arguments}
}
\value{
  A data frame is returned by the function \code{as.data.frame}.

  A single covariate is returned by the function \code{addvar} which
  can be added to an object of class \code{meta}. Internally, the
  \code{merge} function is utilised. See help page
  \code{\link{metagen}} for an example on the use of R function
  addvar.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}, \code{\link{forest.meta}}}

\examples{
data(Fleiss93cont)
#
# Generate additional variable with grouping information
#
Fleiss93cont$group <- c(1,2,1,1,2)
#
# Do meta-analysis without grouping information
#
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, study,
                  data=Fleiss93cont, sm="SMD")
#
# Update meta-analysis object and do subgroup analyses
#
summary(update(meta1, byvar=group))

#
# Same result using metacont function directly
#
meta2 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, study,
                  data=Fleiss93cont, sm="SMD", byvar=group)
summary(meta2)

#
# Compare printout of the following two commands
#
as.data.frame(meta1)
meta1$data
}
\keyword{print}
