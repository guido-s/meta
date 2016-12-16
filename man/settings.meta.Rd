\name{settings.meta}
\alias{settings.meta}

\title{Print and change default settings to conduct and print or plot
meta-analyses in R package \bold{meta}.}

\description{
  Print and change default settings to conduct and print or plot
  meta-analyses in R package \bold{meta}. The following general
  settings are available: \emph{Review Manager 5}, \emph{Journal of
  the American Medical Association}.
}
\usage{
settings.meta(...)
}
\arguments{
  \item{...}{Arguments to change default settings.}
}

\details{
  This function can be used to define defaults for several arguments
  (i.e., assignments using \code{\link{gs}}) of the following R
  functions: \code{\link{metabin}}, \code{\link{metacont}},
  \code{\link{metacor}}, \code{\link{metacr}}, \code{\link{metagen}},
  \code{\link{metainc}}, \code{\link{metaprop}},
  \code{\link{metarate}}

  Furthermore, some of these settings are considered to print
  meta-analysis results using \code{\link{print.meta}} and
  \code{\link{print.summary.meta}}, and to produce forest plots using
  \code{\link{forest.meta}}.

  The function can be used to either change individual settings (see
  Examples) or use one of the following general settings:  
  %%
  \itemize{
   \item \code{settings.meta("revman5")}
   \item \code{settings.meta("jama")}
  }
  The first command can be used to reproduce meta-analyses from
  Cochrane reviews conducted with \emph{Review Manager 5} (RevMan 5,
  \url{http://community.cochrane.org/tools/review-production-tools/revman-5})
  and specifies to use a RevMan 5 layout in forest plots. The second
  command can be used to generate forest plots following instructions
  for authors of the \emph{Journal of the American Medical
  Association}
  (\url{http://jamanetwork.com/journals/jama/pages/instructions-for-authors}).

  RevMan 5 settings, in detail:
  %%
  \tabular{lll}{
   \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
   \code{hakn} \tab FALSE \tab method not available in RevMan 5 \cr
   \code{method.tau} \tab "DL" \tab only available method in RevMan 5 \cr
   \code{tau.common} \tab FALSE \tab common between-study variance in subgroups \cr
   \code{MH.exact} \tab FALSE \tab exact Mantel-Haenszel method \cr
   \code{RR.cochrane} \tab TRUE \tab calculation of risk ratios \cr
   \code{layout} \tab "RevMan5" \tab layout for forest plots \cr
   \code{test.overall} \tab TRUE \tab print information on test of overall effect \cr
   \code{digits.I2} \tab 0 \tab number of digits for I-squared measure \cr
   \code{digits.tau2} \tab 0 \tab number of digits for tau-squared \cr
   \code{CIbracket},  \tab "[" \tab \cr
   \code{CIseparator} \tab ", " \tab print confidence intervals as "\code{[., .]}"
  }
  JAMA settings:
  %%
  \tabular{lll}{
   \bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
   \code{layout} \tab "JAMA" \tab layout for forest plots \cr
   \code{test.overall} \tab TRUE \tab print information on test of overall effect \cr
   \code{digits.I2} \tab 0 \tab number of digits for I-squared measure \cr
   \code{CIbracket},  \tab "(" \tab \cr
   \code{CIseparator} \tab "-" \tab print confidence intervals as "\code{(.-.)}"
  }
  
  A list of all arguments with current settings is printed using the
  command \code{settings.meta("print")}.

  In order to reset all settings of R package \bold{meta} the command
  \code{settings.meta("reset")} can be used.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{gs}}, \code{\link{forest.meta}}}

\examples{
#
# Get listing of current settings
#
settings.meta("print")

#
# Meta-analyses using default settings
#
metabin(10, 20, 15, 20)
metaprop(4, 20)
metabin(10, 20, 15, 20, sm="RD")
metaprop(4, 20, sm="PLN")
#
# Change summary measure for R functions metabin and metaprop
# and store old settings
#
oldset <- settings.meta(smbin="RD", smprop="PLN")
#
metabin(10, 20, 15, 20)
metaprop(4, 20)
#
# Use old settings
#
settings.meta(oldset)

#
# Change level used to calculate confidence intervals
# (99%-CI for studies, 99.9%-CI for pooled effects)
#
metagen(1:3, (2:4)/10, sm="MD")
settings.meta(level=0.99, level.comb=0.999)
metagen(1:3, (2:4)/10, sm="MD")

#
# Always print a prediction interval
#
settings.meta(prediction=TRUE)
metagen(1:3, (2:4)/10, sm="MD")
metagen(4:6, (4:2)/10, sm="MD")

#
# Try to set unknown argument results in a warning
#
settings.meta(unknownarg=TRUE)

#
# Reset to default settings of R package meta
#
settings.meta("reset")
metabin(10, 20, 15, 20)
metaprop(4, 20)
metagen(1:3, (2:4)/10, sm="MD")

#
# Do not back transform results (e.g. print log odds ratios instead of
# odds ratios, print transformed correlations/proportions instead of
# correlations/proportions)
#
settings.meta(backtransf=FALSE)
metabin(10, 20, 15, 20)
metaprop(4, 20)
metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))

#
# Forest plot using RevMan 5 style
#
settings.meta("revman5")
forest(metagen(1:3, (2:4)/10, sm="MD", comb.fixed=FALSE),
       label.left="Favours A", label.right="Favours B",
       colgap.studlab = grid::unit(2, "cm"),
       colgap.forest.left = grid::unit(0.2, "cm"))

#
# Forest plot using JAMA style
#
settings.meta("jama")
forest(metagen(1:3, (2:4)/10, sm="MD", comb.fixed=FALSE),
       label.left="Favours A", label.right="Favours B",
       colgap.studlab = grid::unit(2, "cm"),
       colgap.forest.left = grid::unit(0.2, "cm"))
#
# Use slightly different layout for confidence intervals
# (especially useful if upper confidence limit can be negative)
#
settings.meta(CIseparator=" - ")
forest(metagen(-(1:3), (2:4)/10, sm="MD", comb.fixed=FALSE),
       label.left="Favours A", label.right="Favours B",
       colgap.studlab = grid::unit(2, "cm"),
       colgap.forest.left = grid::unit(0.2, "cm"))

# Use old settings
#
settings.meta(oldset)
}
