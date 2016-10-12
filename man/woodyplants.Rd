\name{woodyplants}
\docType{data}
\alias{woodyplants}
\title{Elevated CO_2 and total biomass of woody plants}

\description{Meta-analysis on effects of elevated CO_2 on total biomass of woody plants}

\usage{data(woodyplants)}

\format{A data frame with the following columns:
  \describe{
  \item{obsno}{Observation number}
  \item{papno}{Database paper number}
  \item{treat}{Treatment code}
  \item{level}{Treatment level}
  \item{n.elev}{Number of observations in experimental group (elevated CO_2-level)}
  \item{mean.elev}{Estimated mean in experimental group}
  \item{sd.elev}{Standard deviation in experimental group}
  \item{n.amb}{Number of observations in control group (ambient CO_2-level)}
  \item{mean.amb}{Estimated mean in control group}
  \item{sd.amb}{Standard deviation in control group}
  }
}

\details{
  This dataset has been used as an example in Hedges et al. (1999) to
  describe methods for the meta-analysis of response ratios. The
  complete dataset with 102 observations and 26 variables is available
  online as a supplement. Here only a subset of 10 variables is
  provided and used in the examples.
}

\source{
  Website \url{http://www.esapubs.org/archive/ecol/E080/008/}
}

\references{
  Hedges LV, Gurevitch J, Curtis PS (1999),
  The meta-analysis of response ratios in experimental ecology.
  \emph{Ecology},
  \bold{80}, 1150--6.
}

\examples{
data(woodyplants)

# Meta-analysis of response ratios (Hedges et al., 1999)
#
meta8 <- metacont(n.elev, mean.elev, sd.elev,
                  n.amb, mean.amb, sd.amb,
                  data=woodyplants, sm="ROM",
                  studlab=paste(obsno, papno, sep = " / "))
summary(meta8, prediction=TRUE)

# Meta-analysis for plants grown with low soil fertility treatment
#
meta9 <- update(meta8, subset=(treat=="fert"&level=="low"))
summary(meta9, prediction=TRUE)

# Meta-analysis for plants grown under low light conditions
#
meta10 <- update(meta8, subset=(treat=="light"&level=="low"))
summary(meta10, prediction=TRUE)
}

\keyword{datasets}
