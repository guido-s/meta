\name{Fleiss93cont}
\docType{data}
\alias{Fleiss93cont}
\title{Mental Health Treatment}

\description{
Meta-analysis on the Effect of Mental Health Treatment on Medical
  Utilisation.

Data example in Fleiss (1993) for meta-analysis with continuous
outcomes.
}

\usage{data(Fleiss93cont)}

\format{A data frame with the following columns:
  \describe{
  \item{study}{Study label}
  \item{year}{Year of publication}
  \item{n.e}{Number of observations in psychotherapy group}
  \item{mean.e}{Estimated mean in psychotherapy group}
  \item{sd.e}{Standard deviation in psychotherapy group}
  \item{n.c}{Number of observations in control group}
  \item{mean.c}{Estimated mean in control group}
  \item{sd.c}{Standard deviation in control group}
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
