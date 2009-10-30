\name{funnel}
\alias{funnel}
\title{Generic function to produce a funnel plot.}
\description{
  Draw a funnel or radial plot to assess funnel plot asymmetry in the
  active graphics window.
  
  A contour-enhanced funnel plot can be produced for assessing causes of
  funnel plot asymmetry.
}
\usage{
funnel(x, y, ...)
}
\arguments{
  \item{x}{An object of class \code{meta}, or estimated treatment
    effect in individual studies.}
  \item{y}{Standard error of estimated treatment effect (mandatory if
    \code{x} not of class \code{meta}).}
  \item{\dots}{Graphical parameters as in \code{par} may also be
    passed as arguments.}
}
\details{
  For simple funnel plots, \code{funnel.default} will be used. For an
  object of class \code{meta} the function \code{funnel.meta} will be
  used instead.
  
  A funnel plot or radial plot, also called Galbraith plot, is drawn in
  the active graphics window. If \code{comb.fixed} is TRUE, the pooled
  estimate of the fixed effect model is plotted. If \code{level} is not
  NULL, the corresponding confidence limits are drawn.

  In the funnel plot, if \code{yaxis} is \code{"se"}, the standard error
  of the treatment estimates is plotted on the y axis which is likely to
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
\seealso{\code{\link{metabias}}, \code{\link{funnel.default}}, \code{\link{funnel.meta}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 studlab=paste(author, year),
                 sm="RR", meth="I")


oldpar <- par(mfrow=c(2, 2))

##
## Funnel plots
##
funnel(meta1)
##
## Same result as code above:
##
funnel(meta1$TE, meta1$seTE, sm="RR")

##
## Funnel plot with confidence intervals,
## fixed effect estimate and contours
##
cc <- funnel(meta1, comb.fixed=TRUE,
             level=0.95, contour=c(0.9, 0.95, 0.99))$col.contour
legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"), fill=cc)
##
## Contour-enhanced funnel plot with user-chosen colours
##
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
