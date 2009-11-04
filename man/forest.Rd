\name{forest}

\alias{forest}

\title{Forest plot (new plot function for objects of class meta)}

\description{
  Draws a forest plot in the active graphics window (using grid
  graphics system).
}

\usage{
forest(x, byvar=x$byvar, bylab=x$bylab,
       print.byvar=x$print.byvar, sortvar, studlab=TRUE,
       level=x$level, level.comb=x$level.comb,
       comb.fixed=x$comb.fixed, comb.random=x$comb.random,
       overall=TRUE,
       text.fixed="Fixed effect model", text.random="Random effects model",
       lty.fixed=2, lty.random=3, xlab=NULL, xlab.pos=ref, xlim,
       allstudies=TRUE,
       weight,
       ref=ifelse(x$sm \%in\% c("RR", "OR", "HR"), 1, 0),
       leftcols=NULL, rightcols=NULL,
       leftlabs=NULL, rightlabs=NULL,
       lab.e=x$label.e, lab.c=x$label.c,
       lab.e.attach.to.col=NULL, lab.c.attach.to.col=NULL,
       lwd=1,
       at=NULL, label=TRUE,
       fontsize=12, boxsize=0.8,
       plotwidth=unit(6, "cm"), colgap=unit(2, "mm"),
       col.i="black", col.by="darkgray")
}

\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{x$TE}). Parameter \code{byvar} can not be
    used if \code{x} is an object of class \code{metacum} or
    \code{metainf}.}
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
  \item{comb.fixed}{A logical indicating whether fixed effect estimate
  should be plotted.}
  \item{comb.random}{A logical indicating whether random effects estimate
  should be plotted.}
  \item{overall}{A logical indicating whether overall summaries should
    be plotted. This parameter is useful in combination with the
    parameter \code{byvar} if summaries should only be plotted on group
    level.}
  \item{text.fixed}{A character string used in the plot to label the
    pooled fixed effect estimate.}
  \item{text.random}{A character string used in the plot to label the
  pooled random effects estimate.}
  \item{lty.fixed}{Line type (pooled fixed effect estimate).}
  \item{lty.random}{Line type (pooled random effects estimate).}
  \item{xlab}{A label for the x axis.}
  \item{xlab.pos}{A numeric specifying the center of the label
    on the x axis.}
  \item{xlim}{The x limits (min,max) of the plot.}
  \item{allstudies}{A logical indicating whether studies with inestimable
    treatment effects should be plotted.}
  \item{weight}{A character string indicating which type of plotting
    symbols is to be used for individual treatment estimates. One of
    missing (see Details), \code{"same"}, \code{"fixed"}, or
    \code{"random"}, can be abbreviated. Plot symbols have the same size
    for all studies or represent study weights from fixed effect or
    random effects model.}
  \item{ref}{A numerical giving the reference value to be plotted as a
    line in the forest plot. No reference line is plotted if parameter
    \code{ref} is equal to \code{NA}.}
  \item{leftcols}{A character vector specifying (additional) columns to
  be plotted on the left side of the forest plot (see Details).}
  \item{rightcols}{A character vector specifying (additional) columns to
  be plotted on the right side of the forest plot (see Details).}
  \item{leftlabs}{A character vector specifying labels for (additional)
  columns on left side of the forest plot (see Details).}
  \item{rightlabs}{A character vector specifying labels for (additional)
  columns on right side of the forest plot (see Details).}
  \item{lab.e}{Label to be used for experimental group in table
    heading.}
  \item{lab.c}{Label to be used for control group in table
    heading.}
  \item{lab.e.attach.to.col}{A character specifying the column name where
    label \code{lab.e} should be attached to in table heading.}
  \item{lab.c.attach.to.col}{A character specifying the column name where
    label \code{lab.c} should be attached to in table heading.}
  \item{lwd}{The line width, see \code{\link{par}}.}
  \item{at}{The points at which tick-marks are to be drawn, see
  \code{grid.xaxis}.}
  \item{label}{A logical value indicating whether to draw the
    labels on the tick marks, or an expression or character vector which
    specify the labels to use. See \code{\link{grid.xaxis}}.}
  \item{fontsize}{The size of text (in points), see \code{\link{gpar}}.}
  \item{boxsize}{A numeric used to increase or decrease the size of
    boxes in the forest plot.}
  \item{plotwidth}{A unit object specifying width of the forest plot.}
  \item{colgap}{A unit object specifying gap between columns printed on
    left and right side of forest plot.}
  \item{col.i}{The colour for individual study results and confidence limits.}
  \item{col.by}{A character specifying colour to print information on
    subgroups.}
}

\details{
  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window. Sub-group analyses are conducted and
  displayed in the plot if \code{byvar} is not missing.

  The \code{forest} function is based on the grid graphics
  system. Therefore, to plot a new figure in an existing graphics
  window, one has to use the \code{\link{grid.newpage}} function. In
  order to print the forest plot, (i) resize the graphics window,
  (ii) either use \code{\link{dev.copy2eps}} or
  \code{\link{dev.copy2pdf}}. For basic forest plots, the
  \code{\link{plot.meta}} function can be used.
  
  Information from object \code{x} is utilised if argument \code{weight}
  is missing. Weights from the fixed effect model are used
  (\code{weight="fixed"}) if parameter \code{x$comb.fixed} is
  \code{TRUE}; weights from the random effects model are used
  (\code{weight="random"}) if parameter \code{x$comb.random}
  is \code{TRUE} and \code{x$comb.fixed} is \code{FALSE}.

  The parameters \code{leftcols} and \code{rightcols} can be used to
  specify columns which are plotted on the left and right side of the
  forest plot, respectively. If these parameters are \code{NULL}, the
  following default columns will be plotted.

  Parameter \code{rightcols}: (i) estimated treatment effect with
  \code{level}-confidence interval, (ii) in addition, weights of the
  fixed and/or random effects model will be given, if
  \code{comb.fixed=TRUE} and/or \code{comb.random=TRUE}. For an object
  of class \code{metacum} or \code{metainf} only the estimated treatment
  effect with \code{level}-confidence interval are plotted.
  
  Parameter \code{leftcols}:
  (i) \code{leftcols=c("studlab", "event.e", "n.e", "event.c", "n.c")}
  for an object of class \code{metabin},
  (ii) \code{leftcols=c("studlab", "n.e", "mean.e", "sd.e", "n.c",
    "mean.c", "sd.c")} for an object of class \code{metacont},
  (iii) \code{leftcols=c("studlab", "TE", "seTE")}
  for an object of class \code{metagen},
  (iv) \code{leftcols=c("studlab", "event", "n")}
  for an object of class \code{metaprop},
  (v) \code{leftcols=c("studlab")}
  for an object of class \code{metacum} or \code{metainf}.

  The parameters \code{leftlabs} and \code{rightlabs} can be used to
  specify column headings which are plotted on left and right side of
  the forest plot, respectively. For certain columns predefined labels
  exist. If the parameters \code{leftlabs} and \code{rightlabs} are
  \code{NULL}, the following default labels will be used: for columns
  \code{c("studlab", "TE", "seTE", "n.e", "n.c", "event.e", "event.c",
  "mean.e", "mean.c", "sd.e", "sd.c", "effect", "ci", "w.fixed",
  "w.random")} the labels \code{c("Study", "TE", "seTE", "Total",
  "Total", "Events", "Events", "Mean", "Mean", "SD", "SD",
  summary measure, level for confidence interval, "W(fixed)",
  "W(random)")}. For additional columns the column name
  will be used as label. It is possible to only provide labels for new
  columns (see Examples).

  If parameters \code{lab.e} and \code{lab.c} are \code{NULL},
  "Experimental" and "Control" are used as labels for experimental and
  control group, respectively.

  For subgroups (argument \code{byvar} not \code{NULL}), results for the
  fixed effect model will be plotted if both arguments \code{comb.fixed}
  and \code{comb.random} are \code{TRUE}. In order to plot results for
  the random effects model within subgroups, use
  \code{comb.fixed==FALSE} and \code{comb.random==TRUE}.

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

\seealso{\code{\link{plot.meta}}, \code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", meth="I",
                 studlab=paste(author, year))


grid.newpage()
##
## Do forest plot
##
forest(meta1, comb.fixed=TRUE, comb.random=TRUE)


grid.newpage()
##
## Change set of columns printed on left side
## of forest plot
##
forest(meta1, comb.fixed=TRUE, comb.random=FALSE,
       leftcols="studlab")


grid.newpage()
##
## 1. Change order of columns on left side
## 2. Attach labels to columns 'event.e' and 'event.c'
##    instead of columns 'n.e' and 'n.c'
##
forest(meta1,
       leftcols=c("studlab", "n.e", "event.e", "n.c", "event.c"),
       lab.e.attach.to.col="event.e",
       lab.c.attach.to.col="event.c",
       comb.fixed=TRUE)


Olkin95$studlab <- paste(Olkin95$author, Olkin95$year)
##
## Add variables 'year' and 'author' to meta-analysis object
##
meta1$year <- addvar(meta1, Olkin95, "year")
meta1$author <- addvar(meta1, Olkin95, "author")

grid.newpage()
##
## Specify column labels only for newly created variables
## 'year' and 'author'
##
forest(meta1,
       leftcols=c("studlab", "event.e", "n.e", "event.c", "n.c",
                  "author", "year"),
       leftlabs=c("Author", "Year of Publ"),
       comb.fixed=TRUE)
}

\keyword{hplot}
