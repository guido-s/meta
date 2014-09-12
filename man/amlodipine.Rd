\name{amlodipine}
\docType{data}
\alias{amlodipine}
\title{Amlodipine for Work Capacity}

\description{Meta-analysis on the effect of amlodipine on
  work capacity.

 This meta-analysis is used as a data example in Hartung and
 Knapp (2001).
}

\usage{data(amlodipine)}

\format{A data frame with the following columns:
  \describe{
  \item{study}{Study label}
  \item{n.amlo}{Number of observations in amlodipine group}
  \item{mean.amlo}{Estimated mean in amlodipine group}
  \item{var.amlo}{Variance in amlodipine group}
  \item{n.plac}{Number of observations in placebo group}
  \item{mean.plac}{Estimated mean in placebo group}
  \item{var.plac}{Variance in placebo group}
  }
}

\source{
  Hartung J & Knapp G (2001),
  On tests of the overall treatment effect in meta-analysis with
  normally distributed responses.
  \emph{Statistics in Medicine}, \bold{20}, 1771--82.
  doi: 10.1002/sim.791 .
}

\seealso{\code{\link{metacont}}}

\examples{
data(amlodipine)
#
m <- metacont(n.amlo, mean.amlo, sqrt(var.amlo),
              n.plac, mean.plac, sqrt(var.plac),
              data=amlodipine, studlab=study)

m.hakn <- metacont(n.amlo, mean.amlo, sqrt(var.amlo),
                   n.plac, mean.plac, sqrt(var.plac),
                   data=amlodipine, studlab=study,
                   hakn=TRUE)

# Results for mean difference - see Table III in Hartung and Knapp (2001) 
#
res.md <- rbind(data.frame(summary(m)$fixed)[c("TE", "lower", "upper")],
                data.frame(summary(m)$random)[c("TE", "lower", "upper")],
                data.frame(summary(m.hakn)$random)[c("TE", "lower", "upper")])
#
res.md <- round(res.md, 5)
#
row.names(res.md) <- c("FE", "RE", "RE (HaKn)")
names(res.md) <- c("Absolute difference", "CI lower", "CI upper")
#
res.md
}

\keyword{datasets}
