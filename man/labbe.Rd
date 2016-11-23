\name{labbe.metabin}

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
      xlab = NULL, ylab = NULL,
      TE.fixed = NULL, TE.random = NULL,
      comb.fixed = !is.null(TE.fixed), comb.random = !is.null(TE.random),
      backtransf = TRUE,
      axes = TRUE,
      pch = 21, text = NULL, cex = 1,
      col = "black", bg = "lightgray",
      lwd = 1, lwd.fixed = lwd, lwd.random = lwd,
      lty.fixed = 2, lty.random = 9,
      col.fixed = col, col.random = col,
      nulleffect = TRUE,
      lwd.nulleffect = lwd, col.nulleffect = "lightgray",
      sm = NULL, weight,
      studlab = FALSE, cex.studlab = 0.8,
      label.e = NULL, label.c = NULL,
      ...)

\method{labbe}{metabin}(x,
      xlim, ylim,
      xlab = NULL, ylab = NULL,
      TE.fixed = x$TE.fixed,
      TE.random = x$TE.random,
      comb.fixed = x$comb.fixed,
      comb.random = x$comb.random,
      backtransf = x$backtransf,
      axes = TRUE,
      pch = 21, text = NULL, cex = 1,
      col = "black", bg = "lightgray",
      lwd = 1, lwd.fixed = lwd, lwd.random = lwd,
      lty.fixed = 2, lty.random = 9,
      col.fixed = col, col.random = col,
      nulleffect = TRUE,
      lwd.nulleffect = lwd, col.nulleffect = "lightgray",
      sm = x$sm, weight,
      studlab = FALSE, cex.studlab = 0.8,
      label.e = x$label.e, label.c = x$label.c,
      ...)
}
\arguments{
  \item{x}{The x coordinates of points of the L'Abbé
    plot. Alternatively, an object of class \code{metabin}.}
  \item{y}{The y coordinates of the L'Abbé plot, optional if \code{x}
    is an appropriate structure.}
\item{xlim}{The x limits (min, max) of the plot.}
  \item{ylim}{The y limits (min, max) of the plot.}
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
  \item{backtransf}{A logical indicating which values should be
    printed on x- and y-axis (see Details).}
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
  \item{lwd.fixed}{The line width(s) for fixed effect estimate(s) (if
    \code{comb.fixed} is not \code{NULL} or \code{FALSE}).}
  \item{lwd.random}{The line width(s) for random effects estimate(s)
    (if \code{comb.random} is not \code{NULL} or \code{FALSE}).}
  \item{lty.fixed}{Line type(s) for fixed effect estimate(s).}
  \item{lty.random}{Line type(s) for random effects estimate(s).}
  \item{col.fixed}{Color of line(s) for fixed effect estimate(s).}
  \item{col.random}{Color of line(s) for random effects estimate(s).}
  \item{nulleffect}{A logical indicating whether line for null effect
    should be added to the plot..}
  \item{lwd.nulleffect}{Width of line for null effect.}
  \item{col.nulleffect}{Color of line for null effect.}
  \item{sm}{A character string indicating underlying summary measure,
    i.e., \code{"RD"}, \code{"RR"}, \code{"OR"}, or \code{"ASD"}.}
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
  \item{label.e}{Label for experimental group.}
  \item{label.c}{Label for control group.}
  \item{\dots}{Graphical arguments as in \code{par} may also be
    passed as arguments.}
}
\details{
  A L'Abbé plot is a scatter plot with the risk in the control group
  on the x-axis and the risk in the experimental group on the y-axis
  (L'Abbé et al., 1987). It can be used to evaluate heterogeneity in
  meta-analysis. Furthermore, this plot can aid to choose a summary
  measure (odds ratio, risk ratio, risk difference) that will result
  in more consistent results.

  If argument \code{backtransf} is TRUE (default), event probabilities
  will be printed on x- and y-axis. Otherwise, transformed event
  probabilities will be printed as defined by the summary measure,
  i.e., log odds of probabilities for odds ratio as summary measure
  (\code{sm = "OR"}), log probabilities for \code{sm = "RR"}, and
  arcsine-transformed probabilities for \code{sm = "ASD"}.
  
  If \code{comb.fixed} is TRUE, the pooled estimate of the fixed effect
  model is plotted as a line. If \code{comb.random} is TRUE, the pooled
  estimate of the random effects model is plotted as a line.

  Information from object \code{x} is utilised if argument
  \code{weight} is missing. Weights from the fixed effect model are
  used (\code{weight = "fixed"}) if argument \code{x$comb.fixed} is
  \code{TRUE}; weights from the random effects model are used
  (\code{weight = "random"}) if argument \code{x$comb.random} is
  \code{TRUE} and \code{x$comb.fixed} is \code{FALSE}.
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
                 data = Olkin95,
                 studlab = paste(author, year),
                 sm = "RR", method = "I")


#
# L'Abbe plot for risk ratio
#
labbe(meta1)

# L'Abbe plot for odds ratio
#
labbe(meta1, sm = "OR")
# same plot
labbe(update(meta1, sm = "OR"))

# L'Abbe plot for risk difference
#
labbe(meta1, sm = "RD")

# L'Abbe plot on log odds scale
#
labbe(meta1, sm = "OR", backtransf = FALSE)

# L'Abbe plot for odds ratio with coloured lines for various treatment
# effects (defined as log odds ratios)
#
mycols <- c("blue", "yellow", "green", "red",
            "green", "yellow", "blue")
labbe(meta1, sm = "OR",
      comb.random = FALSE,
      TE.fixed = log(c(1 / 10, 1 / 5, 1 / 2,
                       1, 2, 5, 10)),
      col.fixed = mycols, lwd.fixed = 2)

# L'Abbe plot on log odds scale with coloured lines for various
# treatment effects (defined as log odds ratios)
#
labbe(meta1, sm = "OR",
      comb.random = FALSE,
      TE.fixed = log(c(1 / 10, 1 / 5, 1 / 2,
                       1, 2, 5, 10)),
      col.fixed = mycols, lwd.fixed = 2,
      backtransf = FALSE)
}

\keyword{hplot}
