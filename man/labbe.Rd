\name{labbe}

\alias{labbe}
\alias{labbe.metabin}
\alias{labbe.default}

\title{L'Abbe plot for meta-analysis with binary outcomes}

\description{
  Draw a L'Abbé plot for meta-analysis with binary outcomes.
}

\usage{
labbe(x, ...)

\method{labbe}{default}(x, y,
      xlim, ylim,
      xlab=NULL, ylab=NULL,
      TE.fixed, TE.random,
      comb.fixed=FALSE, comb.random=FALSE,
      axes=TRUE,
      pch=21, text=NULL, cex=1,
      col="black", bg="lightgray",
      lwd=1, lwd.fixed=lwd, lwd.random=lwd,
      lty.fixed=2, lty.random=9,
      sm=NULL, weight,
      studlab=FALSE, cex.studlab=0.8,
      ...)

\method{labbe}{metabin}(x,
      xlim, ylim,
      xlab=NULL, ylab=NULL,
      TE.fixed=x$TE.fixed,
      TE.random=x$TE.random,
      comb.fixed=x$comb.fixed,
      comb.random=x$comb.random,
      axes=TRUE,
      pch=21, text=NULL, cex=1,
      col="black", bg="lightgray",
      lwd=1, lwd.fixed=lwd, lwd.random=lwd,
      lty.fixed=2, lty.random=9,
      sm=x$sm, weight,
      studlab=FALSE, cex.studlab=0.8, ...)
}
\arguments{
  \item{x}{The x coordinates of points of the L'Abbé
    plot. Alternatively, an object of class \code{metabin}.}
  \item{y}{The y coordinates of the L'Abbé plot, optional if \code{x}
    is an appropriate structure.}
\item{xlim}{The x limits (min,max) of the plot.}
  \item{ylim}{The y limits (min,max) of the plot.}
  \item{xlab}{A label for the x-axis.}
  \item{ylab}{A label for the y-axis.}
  \item{TE.fixed}{A numeric or vector specifying combined fixed effect
    estimate(s).}
  \item{TE.random}{A numeric or vector specifying combined random
    effects estimate(s).}
  \item{comb.fixed}{A logical indicating whether the pooled fixed
    effect estimate should be plotted.}
  \item{comb.random}{A logical indicating whether the pooled random
    effects estimate should be plotted.}
  \item{axes}{A logical indicating whether axes should be drawn on the
    plot.}
  \item{pch}{The plotting symbol used for individual studies.}
  \item{text}{A character vector specifying the text to be used instead
    of plotting symbol.}
  \item{cex}{The magnification to be used for plotting symbol.}
  \item{col}{A vector with colour of plotting symbols.}
  \item{bg}{A vector with background colour of plotting symbols (only
    used if \code{pch} in \code{21:25}).}
  \item{lwd}{The line width.}
  \item{lwd.fixed}{The line width for fixed effect estimate (if
    \code{comb.fixed} is not \code{NULL} or \code{FALSE}).}
  \item{lwd.random}{The line width for random effects estimate (if
    \code{comb.random} is not \code{NULL} or \code{FALSE}).}
  \item{lty.fixed}{Line type (pooled fixed effect estimate).}
  \item{lty.random}{Line type (pooled random effects estimate).}
  \item{sm}{A character string indicating underlying summary measure,
    i.e., \code{"RD"}, \code{"RR"}, \code{"OR"}.}
  \item{weight}{Either a numeric vector specifying relative sizes of
    plotting symbols or a character string indicating which type of
    plotting symbols is to be used for individual treatment
    estimates. One of missing (see Details), \code{"same"},
    \code{"fixed"}, or \code{"random"}, can be abbreviated. Plot symbols
    have the same size for all studies or represent study weights from
    fixed effect or random effects model.}
  \item{studlab}{A logical indicating whether study labels should be
    printed in the graph. A vector with study labels can also be
    provided (must be of same length as \code{x$event.e} then).}
  \item{cex.studlab}{Size of study labels.}
  \item{\dots}{Graphical arguments as in \code{par} may also be
    passed as arguments.}
}
\details{
  A L'Abbé plot is a scatter plot with the risk in the control group
  on the x-axis and the risk in the experimental group on the y-axis
  (L'Abbé et al., 1987).  It can be used to evaluate heterogeneity in
  meta-analysis. Furthermore, this plot can aid to choose a summary
  measure (odds ratio, risk ratio, risk difference) that will result
  in more consistent results.
  
  If \code{comb.fixed} is TRUE, the pooled estimate of the fixed effect
  model is plotted as a line. If \code{comb.random} is TRUE, the pooled
  estimate of the random effects model is plotted as a line.

  Information from object \code{x} is utilised if argument \code{weight}
  is missing. Weights from the fixed effect model are used
  (\code{weight="fixed"}) if argument \code{x$comb.fixed} is
  \code{TRUE}; weights from the random effects model are used
  (\code{weight="random"}) if argument \code{x$comb.random}
  is \code{TRUE} and \code{x$comb.fixed} is \code{FALSE}.
}

\references{
  L'Abbé  KA, Detsky AS, O'Rourke K (1987),
  Meta-analysis in clinical research. \emph{Annals of Internal
  Medicine}, \bold{107}, 224--233.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95,
                 studlab=paste(author, year),
                 sm="RR", method="I")


#
# L'Abbe plot
#
labbe(meta1)
}

\keyword{hplot}
