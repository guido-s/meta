\name{metabind}
\alias{metabind}

\title{Combine meta-analysis objects}

\description{
This function can be used to combine meta-analysis objects and is, for
example, useful to generate a forest plot with results of subgroup
analyses.
}

\usage{
metabind(..., name, pooled, backtransf, outclab)
}

\arguments{
  \item{...}{Any number of meta-analysis objects (see Details).}
  \item{name}{An optional character vector providing descriptive names
    for the meta-analysis objects.}
  \item{pooled}{A character string indicating whether results of a
    fixed effect or random effects model should be considered. Either
    \code{"fixed"} or \code{"random"}, can be abbreviated.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and plots. If \code{backtransf=TRUE}
    (default), results for \code{sm="OR"} are printed as odds ratios
    rather than log odds ratios, for example.}
  \item{outclab}{Outcome label for all meta-analyis objects.}
}

\details{
This function can be used to combine meta-analysis objects and is, for
example, useful to generate a forest plot with results of subgroup
analyses.
}

\value{
  An object of class \code{c("metabind", "meta")} with corresponding
  \code{print}, \code{summary}, and \code{forest} functions. See
  \code{\link{metagen}} for more information on list elements.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metagen}}, \code{\link{forest.metabind}}}

\examples{
data(Fleiss93cont)

# Add some (fictitious) grouping variables:
Fleiss93cont$age <- c(55, 65, 55, 65, 55)
Fleiss93cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")

meta1 <- metacont(n.e, mean.e, sd.e,
                  n.c, mean.c, sd.c,
                  data = Fleiss93cont, sm = "MD")

mu1 <- update(meta1, byvar = age, bylab = "Age group")
mu2 <- update(meta1, byvar = region, bylab = "Region")

mb1 <- metabind(mu1, mu2)
mb1

forest(mb1)
}
