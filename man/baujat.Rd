\name{baujat}

\alias{baujat}
\alias{baujat.meta}

\title{Baujat plot to explore heterogeneity in meta-analysis}

\description{
  Draw a Baujat plot to explore heterogeneity in meta-analysis.
}

\usage{
baujat(x, ...)

\method{baujat}{meta}(x,
       yscale=1, xlim, ylim,
       xlab="Contribution to overall heterogeneity",
       ylab="Influence on overall result",
       pch=21, cex=1, col="black", bg="darkgray",
       studlab=TRUE, cex.studlab=0.8,
       xmin=0, ymin=0, pos=2, offset=0.5,
       grid=TRUE, col.grid="lightgray", lty.grid="dotted", lwd.grid=par("lwd"),
       pty="s", ...)
}
\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{yscale}{Scaling factor for values on y-axis.}
  \item{xlim}{The x limits (min,max) of the plot.}
  \item{ylim}{The y limits (min,max) of the plot.}
  \item{xlab}{A label for the x-axis.}
  \item{ylab}{A label for the y-axis.}
  \item{pch}{The plotting symbol used for individual studies.}
  \item{cex}{The magnification to be used for plotting symbol.}
  \item{col}{A vector with colour of plotting symbols.}
  \item{bg}{A vector with background colour of plotting symbols (only
    used if \code{pch} in \code{21:25}).}
  \item{studlab}{A logical indicating whether study labels should be
    printed in the graph. A vector with study labels can also be
    provided (must be of same length as \code{x$TE} then).}
  \item{cex.studlab}{The magnification for study labels.}
  \item{xmin}{A numeric specifying minimal value to print study labels (on x-axis).}
  \item{ymin}{A numeric specifying minimal value to print study labels (on y-axis).}
  \item{pos}{A position specifier for study labels (see \code{\link{text}}).}
  \item{offset}{Offset for study labels (see \code{\link{text}}).}
  \item{grid}{A logical indicating whether a grid is printed in the plot.}
  \item{col.grid}{Colour for grid lines.}
  \item{lty.grid}{The line type for grid lines.}
  \item{lwd.grid}{The line width for grid lines.}
  \item{pty}{A character specifying type of plot region (see \code{\link{par}}).}
  \item{\dots}{Graphical arguments as in \code{par} may also be
    passed as arguments.}
}
\details{
  Baujat et al. (2002) introduced a scatter plot to explore
  heterogeneity in meta-analysis. On the x-axis the contribution of
  each study to the overall heterogeneity statistic (see list object
  \code{Q} of the meta-analysis object \code{x}) is plotted. On the
  y-axis the standardised difference of the overall treatment effect
  with and without each study is plotted; this quantity describes the
  influence of each study on the overal treatment effect.

  Internally, the \code{\link{metainf}} function is used to calculate
  the values on the y-axis. 
}

\value{
  A data.frame with the following variables:
  \item{x}{Coordinate on x-axis (contribution to heterogeneity statistic).}
  \item{y}{Coordinate on y-axis (influence on overall treatment effect).}
}

\references{
  Baujat B, Mahé C, Pignon JP, Hill C (2002),
  A graphical method for exploring heterogeneity in meta-analyses:
  Application to a meta-analysis of 65 trials.
  \emph{Statistics in Medicine}, \bold{30}, 2641--2652.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metagen}}, \code{\link{metainf}}}

\examples{
data(Olkin95)

m1 <- metabin(event.e, n.e, event.c, n.c, data=Olkin95,
              studlab=author, sm="OR", method="I")

# Generate Baujat plot
baujat(m1)

# Do not print study labels if the x-value is smaller than 4 and the
# y-value is smaller than 1.
baujat(m1, yscale=10, xmin=4, ymin=1)

# Change position of study labels
baujat(m1, yscale=10, xmin=4, ymin=1,
       pos=1, xlim=c(0, 6.5))

# Generate Baujat plot and assign x- and y- coordinates to R object b1
b1 <- baujat(m1)

# Calculate overall heterogeneity statistic
sum(b1$x)
m1$Q
}

\keyword{hplot}
