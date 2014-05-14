\name{Fleiss93}
\docType{data}
\alias{Fleiss93}
\title{Aspirin after Myocardial Infarction}

\description{
Meta-analysis on aspirin in preventing death after myocardial
infarction.

Data example in Fleiss (1993) for meta-analysis with binary outcomes.
}

\usage{data(Fleiss93)}

\format{A data frame with the following columns:
  \describe{
  \item{study}{Study label}
  \item{year}{Year of publication}
  \item{event.e}{Number of deaths in aspirin group}
  \item{n.e}{Number of observations in aspirin group}
  \item{event.c}{Number of deaths in placebo group}
  \item{n.c}{Number of observations in placebo group}
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
