\name{Olkin95}
\docType{data}
\alias{Olkin95}
\title{Thrombolytic Therapy after Acute Myocardial Infarction}

\description{Meta-analysis on Thrombolytic Therapy after Acute Myocardial Infarction}

\usage{data(Olkin95)}

\format{A data frame with the following columns:
  \describe{
  \item{author}{first author}
  \item{year}{year of publication}
  \item{event.e}{number of events in experimental group}
  \item{n.e}{number of observations in experimental group}
  \item{event.c}{number of events in control group}
  \item{n.c}{number of observations in control group}
  }
}

\source{
  Olkin I (1995),
  Statistical and theoretical considerations in meta-analysis.
  \emph{Journal of Clinical Epidemiology}, \bold{48}, 133--146.
}

\examples{
data(Olkin95)
summary(metabin(event.e, n.e, event.c, n.c, data=Olkin95))
}

\keyword{datasets}
