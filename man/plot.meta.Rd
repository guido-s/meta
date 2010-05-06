\name{plot.meta}

\alias{plot.meta}

\title{Plot function for objects of class meta (deprecated)}

\description{
  Draws a forest plot in the active graphics window.
}

\usage{
\method{plot}{meta}(x, byvar=x$byvar, bylab=x$bylab,
     print.byvar=x$print.byvar, 
     sortvar, studlab=TRUE, level=x$level, level.comb=x$level.comb,
     comb.fixed=x$comb.fixed, comb.random=x$comb.random, overall=TRUE,
     text.fixed="Fixed effect model", text.random="Random effects model",
     lty.fixed=2, lty.random=3, xlab=NULL, xlim, ylim, lwd=1, cex=1,
     cex.comb=1.2 * cex, cex.axis=cex, cex.lab=cex,
     log=ifelse(x$sm \%in\% c("RR", "OR", "HR"), "x", ""),
     axes=TRUE, allstudies=TRUE,
     weight=ifelse(comb.random, "random", "fixed"), scale.diamond=1,
     scale.square= 1, col.i="black",
     clim=xlim, arrow.length=0.1,
     ref=ifelse(x$sm \%in\% c("RR", "OR", "HR"), 1, 0),
     ...)
}

\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{x$TE}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
  \item{studlab}{A logical indicating whether study labels should be
    printed in the graph. A vector with study labels can also be
    provided (must be of same length as \code{x$TE} then).}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.fixed}{A logical indicating whether fixed effect
    estimate should be plotted.}
  \item{comb.random}{A logical indicating whether random effects
    estimate should be plotted.}
  \item{overall}{A logical indicating whether overall summaries should
    be plotted. This parameter is useful in combination with the
    parameter \code{byvar} if summaries should only be plotted on group
    level.}
  \item{text.fixed}{A character string used in the plot to label the pooled fixed
    effects estimate.}
  \item{text.random}{A character string used in the plot to label the pooled random
    effects estimate.}
  \item{lty.fixed}{Line type (pooled fixed effect estimate).}
  \item{lty.random}{Line type (pooled random effects estimate).}
  \item{xlab}{A label for the x axis.}
  \item{xlim}{The x limits (min,max) of the plot.}
  \item{ylim}{The y limits (min,max) of the plot.}
  \item{lwd}{The line width.}
  \item{cex}{A numerical value giving the amount by which plotting text
    and symbols should be scaled relative to the default.}
  \item{cex.comb}{A numerical value giving the amount by which plotting
    text and symbols for pooled fixed and random effects estimates should be scaled.}
  \item{cex.axis}{The magnification to be used for axis annotation
    relative to the current setting of \code{cex}.}
  \item{cex.lab}{The magnification to be used for x and y labels
    relative to the current setting of \code{cex}.}
  \item{log}{A character string which contains \code{"x"} if the x axis
    is to be logarithmic (other values for \code{log} are not
    reasonable).}
  \item{axes}{A logical indicating whether the x axis should be drawn on
    the plot.}
  \item{allstudies}{A logical indicating whether studies with inestimable
    treatment effects should be plotted.}
  \item{weight}{A character string indicating which type of plotting
    symbols is to be used for individual treatment estimates. One of
    \code{"same"}, \code{"fixed"}, or \code{"random"}, can be
    abbreviated. Plot symbols have the same size for all studies or
    represent study weights from fixed effect or random
    effects model.}
  \item{scale.diamond}{A numerical value giving the amount by which the
    diamond representing pooled treatment effects should be scaled
    relative to the default.}
  \item{scale.square}{A numerical value giving the amount by which the
    square representing treatment effects in individual studies should
    be scaled relative to the default.}
  \item{\dots}{Graphical parameters as in \code{par} may also be
    passed as arguments.}
  \item{col.i}{The colour for individual study results and confidence limits.}
  \item{clim}{Limits (min,max) where to cut confidence limits; arrows
    are plotted if confidence limits are outside the range of clim.}
  \item{arrow.length}{Length of the edges of the arrow head (in inches)
    which is plotted if confidence limits are outside the range of
    clim. See also function \code{arrows}.}
  \item{ref}{A numerical value defining a reference value which is
    plotted as a vertical line.}
}

\details{
  Please note, the R code of the function \code{plot.meta} is no longer
  maintained. Use of the function might result in (unexpected) error
  messages or even in wrong results. Please use the function
  \code{\link{forest.meta}} instead.

  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window. Sub-group analyses are conducted and displayed
  in the plot if \code{byvar} is not missing.
  
  Review Manager 5 (RevMan 5) is the current software used for
  preparing and maintaining Cochrane Reviews
  (\url{http://www.cc-ims.net/revman/}). In RevMan 5, subgroup analyses
  can be defined and data from a Cochrane review can be imported to R using the
  function \code{read.rm5}. If a meta-analysis is then conducted using
  function \code{metacr}, information on subgroups is available in R
  (components \code{byvar}, \code{bylab}, and \code{print.byvar},
  \code{byvar} in an object of class \code{"meta"}). Accordingly, by using function
  \code{metacr} there is no need to define subgroups in order to redo
  the statistical analysis conducted in the Cochrane review.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{forest.meta}}, \code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", meth="I")

oldpar <- par(mfrow=c(2, 2))

plot(meta1)
plot(meta1, byvar=c(1,2,1,2), bylab="label")
plot(meta1, byvar=1:4, xlim=c(0.02, 10))

par(oldpar)
}

\keyword{hplot}
