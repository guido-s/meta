\name{bubble}

\alias{bubble}
\alias{bubble.metareg}

\title{Bubble plot to display the result of a meta-regression}

\description{
  Draw a bubble plot to display the result of a meta-regression.
}

\usage{
bubble(x, ...)

\method{bubble}{metareg}(x,
       xlim, ylim, xlab, ylab,
       cex, min.cex=0.5, max.cex=5,
       pch=21, col="black", bg="darkgray",
       lty=1, lwd=1, col.line="black",
       studlab=FALSE, cex.studlab=0.8,
       pos=2, offset=0.5,
       regline=TRUE,
       axes=TRUE, box=TRUE,
       ...)
}
\arguments{
  \item{x}{An object of class \code{metareg}.}
  \item{xlim}{The x limits (min,max) of the plot.}
  \item{ylim}{The y limits (min,max) of the plot.}
  \item{xlab}{A label for the x-axis.}
  \item{ylab}{A label for the y-axis.}
  \item{cex}{The magnification to be used for plotting symbols.}
  \item{min.cex}{Minimal magnification for plotting symbols.}
  \item{max.cex}{Maximal magnification for plotting symbols.}
  \item{pch}{The plotting symbol used for individual studies.}
  \item{col}{A vector with colour of plotting symbols.}
  \item{bg}{A vector with background colour of plotting symbols (only
    used if \code{pch} in \code{21:25}).}
  \item{lty}{The line type for the meta-regression line.}
  \item{lwd}{The line width for the meta-regression line.}
  \item{col.line}{Colour for the meta-regression line.}
  \item{studlab}{A logical indicating whether study labels should be
    printed in the graph. A vector with study labels can also be
    provided (must be of same length as the numer of studies in the
    meta-analysis then).}
  \item{cex.studlab}{The magnification for study labels.}
  \item{pos}{A position specifier for study labels (see \code{\link{text}}).}
  \item{offset}{Offset for study labels (see \code{\link{text}}).}
  \item{regline}{A logical indicating whether a regression line should
    be added to the bubble plot.}
  \item{axes}{A logical indicating whether axes should be printed.}
  \item{box}{A logical indicating whether a box should be printed.}
  \item{\dots}{Graphical arguments as in \code{par} may also be
    passed as arguments.}
}

\details{
  A bubble plot can be used to display the result of a
  meta-regression. It is a scatter plot with the treatment effect for
  each study on the y-axis and the covariate used in the
  meta-regression on the x-axis. Typically, the size of the plotting
  symbol is inversely proportional to the variance of the estimated
  treatment effect (Thompson & Higgins, 2002).

  Argument \code{cex} specifies the plotting size for each individual
  study. If this argument is missing the weights from the
  meta-regression model will be used (which typically is a random
  effects model). Use \code{weight="fixed"} in order to utilise
  weights from a fixed effect model to define the size of the plotted
  symbols (even for a random effects meta-regression). If a vector
  with individual study weights is provided, the length of this vector
  must be of the same length as the number of studies.

  Arguments \code{min.cex} and \code{max.cex} can be used to define
  the size of the smallest and largest plotting symbol. The plotting
  size of the most precise study is set to \code{max.cex} whereas the
  plotting size of all studies with a plotting size smaller than
  \code{min.cex} will be set to \code{min.cex}.

  For a meta-regression with more than one covariate. Only a scatter
  plot of the first covariate in the regression model is shown. In
  this case the effect of the first covariate adjusted for other
  covariates in the meta-regression model is shown.

  For a factor or categorial covariate separate bubble plots for each
  group compared to the baseline group are plotted.
}

\references{
  Thompson SG, Higgins JP (2002),
  How should meta-regression analyses be undertaken and interpreted?
  \emph{Statistics in Medicine}, \bold{21}, 1559--1573.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metagen}}, \code{\link{metainf}}}

\examples{
data(Fleiss93cont)

# Add some (fictious) grouping variables:
Fleiss93cont$age <- c(55, 65, 52, 65, 58)
Fleiss93cont$region <- c("Europe", "Europe", "Asia", "Asia", "Europe")

meta1 <- metacont(n.e, mean.e, sd.e,
                  n.c, mean.c, sd.c,
                  data=Fleiss93cont, sm="MD")

mr1 <- metareg(meta1, region)
mr1

bubble(mr1)
bubble(mr1, lwd=2, col.line="blue")

mr2 <- metareg(meta1, age)
mr2

bubble(mr2, lwd=2, col.line="blue", xlim=c(50, 70))
bubble(mr2, lwd=2, col.line="blue", xlim=c(50, 70), cex="fixed")

# Do not print regression line
#
bubble(mr2, lwd=2, col.line="blue", xlim=c(50, 70), regline=FALSE)
}

\keyword{hplot}
