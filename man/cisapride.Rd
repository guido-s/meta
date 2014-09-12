\name{cisapride}
\docType{data}
\alias{cisapride}
\title{Cisapride in Non-Ulcer Dispepsia}

\description{Meta-analysis on cisapride in non-ulcer dispepsia.

 This meta-analysis is used as a data example in Hartung and
 Knapp (2001).
}

\usage{data(cisapride)}

\format{A data frame with the following columns:
  \describe{
  \item{study}{Study label}
  \item{event.cisa}{Number of events in cisapride group}
  \item{n.cisa}{Number of observations in cisapride group}
  \item{event.plac}{Number of events in placebo group}
  \item{n.plac}{Number of observations in placebo group}
  }
}

\source{
  Hartung J & Knapp G (2001),
  A Refined Method for the Meta-analysis of Controlled Clinical Trials
  with Binary Outcome.
  \emph{Statistics in Medicine}, \bold{20}, 3875--89.
}

\seealso{\code{\link{metabin}}}

\examples{
data(cisapride)

m.or <- metabin(event.cisa, n.cisa, event.plac, n.plac,
                data=cisapride, sm="OR", method="Inverse",
                studlab=study, addincr=TRUE)

m.rr <- metabin(event.cisa, n.cisa, event.plac, n.plac,
                data=cisapride, sm="RR", method="Inverse",
                studlab=study, addincr=TRUE)

m.or.hakn <- metabin(event.cisa, n.cisa, event.plac, n.plac,
                     data=cisapride, sm="OR", method="Inverse",
                     studlab=study, addincr=TRUE,
                     hakn=TRUE)

m.rr.hakn <- metabin(event.cisa, n.cisa, event.plac, n.plac,
                     data=cisapride, sm="RR", method="Inverse",
                     studlab=study, addincr=TRUE,
                     hakn=TRUE)

# Results for log risk ratio - see Table VII in Hartung and Knapp (2001) 
#
res.rr <- rbind(data.frame(summary(m.rr)$fixed)[c("TE", "lower", "upper")],
                data.frame(summary(m.rr)$random)[c("TE", "lower", "upper")],
                data.frame(summary(m.rr.hakn)$random)[c("TE", "lower", "upper")])
#
row.names(res.rr) <- c("FE", "RE", "RE (HaKn)")
names(res.rr) <- c("Log risk ratio", "CI lower", "CI upper")
#
res.rr


# Results for log odds ratio (Table VII in Hartung and Knapp 2001) 
#
res.or <- rbind(data.frame(summary(m.or)$fixed)[c("TE", "lower", "upper")],
                data.frame(summary(m.or)$random)[c("TE", "lower", "upper")],
                data.frame(summary(m.or.hakn)$random)[c("TE", "lower", "upper")])
#
row.names(res.or) <- c("FE", "RE", "RE (HaKn)")
names(res.or) <- c("Log odds ratio", "CI lower", "CI upper")
#
res.or
}

\keyword{datasets}
