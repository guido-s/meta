\name{Fleiss93cont}
\docType{data}
\alias{Fleiss93cont}
\title{Mental Health Treatment}

\description{Meta-analysis on the Effect of Mental Health Treatment on
  Medical Utilisation}

\usage{data(Fleiss93cont)}

\format{A data frame with the following columns:
  \describe{
  \item{study}{study label}
  \item{year}{year of publication}
  \item{n.e}{number of observations in experimental group}
  \item{mean.e}{estimated mean in experimental group}
  \item{sd.e}{standard deviation in experimental group}
  \item{n.c}{number of observations in control group}
  \item{mean.c}{estimated mean in control group}
  \item{sd.c}{standard deviation in control group}
  }
}

\source{
  Fleiss JL (1993),
  The statistical basis of meta-analysis.
  \emph{Statistical Methods in Medical Research}, \bold{2}, 121--145.
}

\seealso{\code{\link{Fleiss93}}}

\examples{
data(Fleiss93cont)
metacont(n.e, mean.e, sd.e,
         n.c, mean.c, sd.c,
         data=Fleiss93cont,
         studlab=paste(study, year),
         comb.random=FALSE)
}

\keyword{datasets}
