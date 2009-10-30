\name{Fleiss93}
\docType{data}
\alias{Fleiss93}
\title{Aspirin after Myocardial Infarction}

\description{Meta-analysis on Aspirin in Preventing Death after Myocardial Infarction}

\usage{data(Fleiss93)}

\format{A data frame with the following columns:
  \describe{
  \item{study}{study label}
  \item{year}{year of publication}
  \item{event.e}{number of events in experimental group}
  \item{n.e}{number of observations in experimental group}
  \item{event.c}{number of events in control group}
  \item{n.c}{number of observations in control group}
  }
}

\source{
  Fleiss JL (1993),
  The statistical basis of meta-analysis.
  \emph{Statistical Methods in Medical Research}, \bold{2}, 121--145.
}

\examples{
data(Fleiss93)
metabin(event.e, n.e, event.c, n.c,
        data=Fleiss93,
        studlab=paste(study, year),
        sm="OR", comb.random=FALSE)
}

\keyword{datasets}
