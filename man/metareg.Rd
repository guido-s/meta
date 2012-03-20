\name{metareg}

\alias{metareg}

\title{Meta-regression}

\description{
  Meta-regression for objects of class \code{meta}. This is a wrapper
  function for the R function \code{rma.uni} in the R package metafor
  (Viechtbauer 2010).
}

\usage{
metareg(formula, data, method.tau=data$method.tau)
}

\arguments{
  \item{formula}{A formula object.}
  \item{data}{An object of class \code{meta}.}
  \item{method.tau}{A character string indicating which method is used
    to estimate the between-study variance tau-squared. Either
    \code{"DL"}, \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
    \code{"HE"}, or \code{"EB"}, can be abbreviated.}
}

\details{
  This R function is a wrapper function for R function \code{rma.uni}
  in the R package metafor (Viechtbauer 2010), i.e. the function
  \code{metareg} can only be used if the R package metafor is
  installed.
}

\value{
An object of class \code{c("rma.uni","rma")}. Please look at the help page of R function \code{rma.uni} for more details.}

\references{
  Viechtbauer W (2010),
  Conducting Meta-Analyses in R with the Metafor Package.
  \emph{Journal of Statistical Software}, \bold{36}, 1--48.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{summary.meta}}, \code{\link{metagen}}}

\examples{
data(Fleiss93cont)
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont, sm="MD")
meta1$grp <- c(1,2,1,1,2)
metareg(~grp, data=meta1)
}

\keyword{print}
