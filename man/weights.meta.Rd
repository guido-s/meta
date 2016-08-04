\name{weights.meta}

\alias{weights.meta}

\title{Calculate absolute and percentage weights for meta-analysis}

\description{
  The \code{weights.meta} method returns a data frame containing
  information on absolute and percentage weights of individual studies
  contributing to fixed effect and random effects meta-analysis.
}

\usage{
\method{weights}{meta}(object, comb.fixed = object$comb.fixed,
        comb.random = object$comb.random, ...)
}

\arguments{
  \item{object}{An object of class \code{meta}.}
  \item{comb.fixed}{A logical indicating whether absolute and
    percentage weights from the fixed effect model should be
    calculated.}
  \item{comb.random}{A logical indicating whether absolute and
    percentage weights from the random effects model should be
    calculated.}
  \item{\dots}{Additional arguments (ignored at the moment).}
}

\value{
  A data frame with the following variables is returned:
  \tabular{lll}{
   \bold{Variable} \tab \bold{Definition} \tab \bold{Condition} \cr
   w.fixed \tab absolute weights in fixed effect model \tab (if \code{comb.fixed = TRUE}) \cr
   p.fixed \tab percentage weights in fixed effect model \tab (if \code{comb.fixed = TRUE}) \cr
   w.random \tab absolute weights in random effects model \tab (if \code{comb.random = TRUE}) \cr
   p.random \tab percentage weights in random effects model \tab (if \code{comb.random = TRUE}) \cr
  }
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
data(Fleiss93cont)
#
# Do meta-analysis (fixed effect and random effects model)
#
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, study,
                  data=Fleiss93cont, sm="SMD")
#
# Print weights for fixed effect and random effects meta-analysis
#
weights(meta1)

#
# Do meta-analysis (only random effects model)
#
meta2 <- update(meta1, comb.fixed = FALSE)
#
# Print weights for random effects meta-analysis
#
weights(meta2)
#
# Print weights for fixed effect and random effects meta-analysis
#
weights(meta2, comb.fixed = TRUE)
}
