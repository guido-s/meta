\name{funnel}

\alias{funnel}
\alias{funnel.meta}
\alias{funnel.default}
\alias{radial}
\alias{galbraith}
\alias{radial.meta}
\alias{radial.default}

\title{Plot to assess funnel plot asymmetry}

\description{
  Draw a funnel plot or radial plot (also called Galbraith plot) to
  assess funnel plot asymmetry in the active graphics window.
  
  A contour-enhanced funnel plot can be produced for assessing causes of
  funnel plot asymmetry.
}

\usage{
funnel(x, ...)

radial(x, ...)

\method{funnel}{default}(x, y,
       xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
       comb.fixed=FALSE, comb.random=FALSE,
       axes=TRUE,
       pch=21, text=NULL, cex=1,
       lty.fixed=2, lty.random=9,
       lwd=1, lwd.fixed=lwd, lwd.random=lwd,
       col="black", bg="darkgray",
       col.fixed="black", col.random="black",
       log="", yaxis="se", sm="",
       contour.levels=NULL, col.contour,
       ref=ifelse(backtransf & is.relative.effect(sm), 1, 0),
       level=NULL,
       studlab=FALSE, cex.studlab=0.8, backtransf=TRUE, ...)

\method{funnel}{meta}(x,
       xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
       comb.fixed=x$comb.fixed, comb.random=x$comb.random,
       axes=TRUE,
       pch=if (!inherits(x, "trimfill")) 21 else ifelse(x$trimfill, 1, 21),
       text=NULL, cex=1,
       lty.fixed=2, lty.random=9,
       lwd=1, lwd.fixed=lwd, lwd.random=lwd,
       col="black", bg="darkgray",
       col.fixed="black", col.random="black",
       log="", yaxis="se",
       contour.levels=NULL, col.contour,
       ref=ifelse(backtransf & is.relative.effect(x$sm), 1, 0),
       level=x$level,
       studlab=FALSE, cex.studlab=0.8, backtransf=x$backtransf, ...)

\method{radial}{default}(x, y, xlim=NULL, ylim=NULL,
       xlab="Inverse of standard error",
       ylab="Standardised treatment effect (z-score)",
       comb.fixed=TRUE, axes=TRUE,
       pch=1, text=NULL, cex=1, col=NULL,
       level=NULL, ...)

\method{radial}{meta}(x, xlim=NULL, ylim=NULL,
       xlab="Inverse of standard error",
       ylab="Standardised treatment effect (z-score)",
       comb.fixed=TRUE, axes=TRUE,
       pch=1, text=NULL, cex=1, col=NULL,
       level=NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{meta}, or estimated treatment
    effect in individual studies.}
  \item{y}{Standard error of estimated treatment effect.}
  \item{xlim}{The x limits (min,max) of the plot.}
  \item{ylim}{The y limits (min,max) of the plot.}
  \item{xlab}{A label for the x-axis.}
  \item{ylab}{A label for the y-axis.}
  \item{comb.fixed}{A logical indicating whether the pooled fixed effect
    estimate should be plotted.}
  \item{comb.random}{A logical indicating whether the pooled random
    effects estimate should be plotted.}
  \item{axes}{A logical indicating whether axes should be drawn on the
    plot.}
  \item{pch}{The plotting symbol used for individual studies.}
  \item{text}{A character vector specifying the text to be used instead
    of plotting symbol.}
  \item{cex}{The magnification to be used for plotting symbol.}
  \item{lty.fixed}{Line type (pooled fixed effect estimate).}
  \item{lty.random}{Line type (pooled random effects estimate).}
  \item{col}{A vector with colour of plotting symbols.}
  \item{bg}{A vector with background colour of plotting symbols (only
    used if \code{pch} in \code{21:25}).}
  \item{col.fixed}{Color of line representign fixed effect estimate.}
  \item{col.random}{Color of line representign random effects estimate.}
  \item{lwd}{The line width for confidence intervals (if \code{level} is
    not \code{NULL}).}
  \item{lwd.fixed}{The line width for fixed effect estimate (if
    \code{comb.fixed} is not \code{NULL}).}
  \item{lwd.random}{The line width for random effects estimate (if
    \code{comb.random} is not \code{NULL}).}
  \item{log}{A character string which contains \code{"x"} if the x-axis
    is to be logarithmic, \code{"y"} if the y-axis is to be logarithmic
    and \code{"xy"} or \code{"yx"} if both axes are to be logarithmic
    (applies only to function \code{funnel}).}
  \item{yaxis}{A character string indicating which type of weights are
    to be used. Either \code{"se"}, \code{"invvar"}, \code{"invse"}, or \code{"size"}
    (applies only to function \code{funnel}). }
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"}, \code{"MD"},
    \code{"SMD"} (applies only to function \code{funnel}).}
  \item{contour.levels}{A numeric vector specifying contour levels to
    produce contour-enhanced funnel plot.}
  \item{col.contour}{Colour of contours.}
  \item{ref}{Reference value (null effect) used to produce
    contour-enhanced funnel plot.}
  \item{level}{The confidence level utilised in the plot. For the funnel
  plot, confidence limits are not drawn if \code{yaxis="size"}.}
  \item{studlab}{A logical indicating whether study labels should be
    printed in the graph. A vector with study labels can also be
    provided (must be of same length as \code{x$TE} then).}
  \item{cex.studlab}{Size of study labels.}
  \item{backtransf}{A logical indicating whether results for relative
    summary measures (argument \code{sm} equal to \code{"OR"},
    \code{"RR"}, \code{"HR"}, or \code{"IRR"}) should be back
    transformed in funnel plots. If \code{backtransf=TRUE}, results
    for \code{sm="OR"} are printed as odds ratios rather than log odds
    ratios, for example.}
  \item{\dots}{Graphical arguments as in \code{par} may also be
    passed as arguments.}
}
\details{
  A funnel plot or radial plot, also called Galbraith plot, is drawn in
  the active graphics window. If \code{comb.fixed} is TRUE, the pooled
  estimate of the fixed effect model is plotted. If \code{level} is not
  NULL, the corresponding confidence limits are drawn.
  
  In the funnel plot, if \code{yaxis} is \code{"se"}, the standard error
  of the treatment estimates is plotted on the y-axis which is likely to
  be the best choice (Sterne & Egger, 2001). Other possible choices for
  \code{yaxis} are \code{"invvar"} (inverse of the variance),
  \code{"invse"} (inverse of the standard error), and \code{"size"}
  (study size).
  
  For \code{yaxis!="size"}, contour-enhanced funnel plots can be
  produced (Peters et al., 2008) by specifying the contour levels
  (argument \code{contour.levels}). By default (argument
  \code{col.contour} missing), suitable gray levels will be used to
  distinguish the contours. Different colours can be chosen by argument
  \code{col.contour}.
  
}
\references{
  Galbraith RF (1988a),
  Graphical display of estimates having differing standard
  errors. \emph{Technometrics}, \bold{30}, 271--281.

  Galbraith RF (1988b),
  A note on graphical presentation of estimated odds ratios from several
  clinical trials. \emph{Statistics in Medicine}, \bold{7}, 889--894.
  
  Light RJ & Pillemer DB (1984),
  \emph{Summing Up. The Science of Reviewing Research}.
  Cambridge: Harvard University Press.

  Peters JL, Sutton AJ, Jones DR, Abrams KR, Rushton L (2008),
  Contour-enhanced meta-analysis funnel plots help distinguish
  publication bias from other causes of asymmetry. \emph{Journal of
    Clinical Epidemiology}, \bold{61}, 991--996.

  Sterne JAC & Egger M (2001),
  Funnel plots for detecting bias in meta-analysis: Guidelines on choice
  of axis. \emph{Journal of Clinical Epidemiology}, \bold{54},
  1046--1055.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}, Petra Graham \email{pgraham@efs.mq.edu.au}}
\seealso{\code{\link{metabias}}, \code{\link{metabin}}, \code{\link{metagen}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 studlab=paste(author, year),
                 sm="RR", method="I")


#
# Radial plot
#
radial(meta1, level=0.95)


oldpar <- par(mfrow=c(2, 2))

#
# Funnel plots
#
funnel(meta1)
#
# Same result as code above:
#
funnel(meta1$TE, meta1$seTE, sm="RR",
       comb.fixed=TRUE, level=0.95)
#
# Funnel plot with confidence intervals,
# fixed effect estimate and contours
#
cc <- funnel(meta1, comb.fixed=TRUE,
             level=0.95, contour=c(0.9, 0.95, 0.99))$col.contour
legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"), fill=cc)
#
# Contour-enhanced funnel plot with user-chosen colours
#
funnel(meta1, comb.fixed=TRUE,
       level=0.95, contour=c(0.9, 0.95, 0.99),
       col.contour=c("darkgreen", "green", "lightgreen"),
       lwd=2, cex=2, pch=16, studlab=TRUE, cex.studlab=1.25)
legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"),
       fill=c("darkgreen", "green", "lightgreen"))

par(oldpar)
}
\keyword{hplot}
