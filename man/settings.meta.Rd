% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.meta.R
\name{settings.meta}
\alias{settings.meta}
\title{Print and change default meta-analysis settings in R package \bold{meta}}
\usage{
settings.meta(..., quietly = TRUE)
}
\arguments{
\item{...}{Arguments to change default settings.}

\item{quietly}{A logical indicating whether information on settings
should be printed.}
}
\description{
Print and change default settings to conduct and print or plot
meta-analyses in R package \bold{meta}. The following general
settings are available: \emph{Review Manager 5}, \emph{Journal of
the American Medical Association}.
}
\details{
This function can be used to define defaults for several arguments
(i.e., assignments using \code{\link{gs}}) of the following R
functions: \code{\link{metabin}}, \code{\link{metacont}},
\code{\link{metacor}}, \code{\link{metacr}}, \code{\link{metagen}},
\code{\link{metainc}}, \code{\link{metaprop}},
\code{\link{metarate}}

Furthermore, some of these settings are considered to print
meta-analysis results and to produce forest plots.

The function can be used to either change individual settings (see
Examples) or use one of the following general settings:
\itemize{
\item \code{settings.meta("RevMan5")}
\item \code{settings.meta("BMJ")}
\item \code{settings.meta("JAMA")}
\item \code{settings.meta("IQWiG5")}
\item \code{settings.meta("IQWiG6")}
\item \code{settings.meta("geneexpr")}
\item \code{settings.meta("IVhet")}
\item \code{settings.meta("meta4")}
\item \code{settings.meta("meta7")}
}

The first command can be used to reproduce meta-analyses from
Cochrane reviews conducted with \emph{Review Manager 5} (RevMan 5)
and specifies to use a RevMan 5 layout in forest plots.

The second command can be used to generate forest plots in BMJ layout.

The third command can be used to generate forest plots following
instructions for authors of the \emph{Journal of the American
Medical Association}. Study labels according to JAMA guidelines can be
generated using \code{\link{labels.meta}}.

The next commands implement the recommendations of the Institute
for Quality and Efficiency in Health Care, Germany (IQWiG)
accordinging to General Methods 5 and 6, respectively
(\url{https://www.iqwig.de/en/about-us/methods/methods-paper/}).

The setting \code{"geneexpr"} can be used to print p-values in
scientific notation and to suppress the calculation of confidence
intervals for the between-study variance.

The setting \code{"IVhet"} can be used for the inverse variance
heterogeneity model (Doi et al., 2015).

The last settings use the default settings of R package
\bold{meta}, version 4 and 7.0-0, respectively, or below.

RevMan 5 settings, in detail:
\tabular{lll}{
\bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
\code{method.random.ci} \tab "classic" \tab only available method
  in RevMan 5 \cr
\code{method.tau} \tab "DL" \tab only available method in RevMan 5
  \cr
\code{method.I2} \tab "Q" \tab only available method in RevMan 5 \cr
\code{tau.common} \tab FALSE \tab common between-study variance in
  subgroups \cr
\code{MH.exact} \tab FALSE \tab exact Mantel-Haenszel method \cr
\code{RR.Cochrane} \tab TRUE \tab calculation of risk ratios \cr
\code{Q.Cochrane} \tab TRUE \tab calculation of heterogeneity
  statistic \cr
\code{exact.smd} \tab FALSE \tab exact formulae for Hedges' g and
  Cohen's d \cr
\code{layout} \tab "RevMan5" \tab layout for forest plots \cr
\code{prediction} \tab FALSE \tab no prediction interval \cr
\code{test.overall} \tab TRUE \tab print information on test of
  overall effect \cr
\code{test.subgroup} \tab TRUE \tab print information on test for
  subgroup differences \cr
\code{test.effect.subgroup} \tab TRUE \tab print information on
  test for effect in subgroups \cr
\code{forest.I2} \tab TRUE \tab show heterogeneity statistic I2 in
  forest plots \cr
\code{forest.tau2} \tab TRUE \tab show between-study heterogeneity \cr
  \tab \tab variance in forest plots \cr
\code{forest.tau} \tab FALSE \tab do not show between-study heterogeneity \cr
  \tab \tab standard deviation in forest plots \cr
\code{forest.Q} \tab TRUE \tab show heterogeneity statistic Q in
  forest plots
  \cr
\code{forest.pval.Q} \tab TRUE \tab show p-value of test for heterogeneity
  in forest plots
  \cr
\code{forest.Rb} \tab FALSE \tab do not show heterogeneity statistic Rb in
  forest plots
  \cr
\code{digits.tau2} \tab 3 \tab number of digits for tau-squared \cr
\code{digits.tau} \tab 4 \tab number of digits for square root of
  tau-squared \cr
\code{digits.I2} \tab 0 \tab number of digits for I-squared measure
  \cr
\code{CIbracket}, \tab "[" \tab \cr
\code{CIseparator} \tab ", " \tab print confidence intervals as
  "\code{[., .]}" \cr
\code{header.line}, \tab TRUE \tab print header line
}

BMJ settings:
\tabular{lll}{
\bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
\code{layout} \tab "BMJ" \tab layout for forest plots \cr
\code{test.overall} \tab TRUE \tab print information on test of
  overall effect \cr
\code{test.subgroup} \tab FALSE \tab print information on test for
  subgroup differences \cr
\code{test.effect.subgroup} \tab FALSE \tab print information on
  test for effect in subgroups \cr
\code{forest.I2} \tab TRUE \tab show heterogeneity statistic I2 in
  forest plots
  \cr
\code{forest.tau2} \tab TRUE \tab show between-study heterogeneity \cr
  \tab \tab variance in forest plots \cr
\code{forest.tau} \tab FALSE \tab do not show between-study heterogeneity \cr
  \tab \tab standard deviation in forest plots \cr
\code{forest.Q} \tab TRUE \tab show heterogeneity statistic Q in
  forest plots
  \cr
\code{forest.pval.Q} \tab TRUE \tab show p-value of test for heterogeneity
  in forest plots
  \cr
\code{forest.Rb} \tab FALSE \tab do not show heterogeneity statistic Rb in
  forest plots
  \cr
\code{digits.I2} \tab 0 \tab number of digits for I-squared measure
  \cr
\code{digits.pval} \tab 2 \tab number of digits for p-values \cr
\code{CIbracket}, \tab "(" \tab \cr
\code{CIseparator} \tab " to " \tab print confidence intervals as
  "\code{(. to .)}" \cr
\code{hetlab}, \tab \tab "Test for heterogeneity: " \cr
\code{header.line}, \tab TRUE \tab print header line
}

JAMA settings:
\tabular{lll}{
\bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
\code{layout} \tab "JAMA" \tab layout for forest plots \cr
\code{test.overall} \tab TRUE \tab print information on test of
  overall effect \cr
\code{test.subgroup} \tab FALSE \tab print information on test for
  subgroup differences \cr
\code{test.effect.subgroup} \tab FALSE \tab print information on
  test for effect in subgroups \cr
\code{forest.I2} \tab TRUE \tab show heterogeneity statistic I2 in
  forest plots
  \cr
\code{forest.tau2} \tab FALSE \tab do not show between-study heterogeneity \cr
  \tab \tab variance in forest plots \cr
\code{forest.tau} \tab FALSE \tab do not show between-study heterogeneity \cr
  \tab \tab standard deviation in forest plots \cr
\code{forest.Q} \tab TRUE \tab show heterogeneity statistic Q in
  forest plots
  \cr
\code{forest.pval.Q} \tab TRUE \tab show p-value of test for heterogeneity
  in forest plots
  \cr
\code{forest.Rb} \tab FALSE \tab do not show heterogeneity statistic Rb in
  forest plots
  \cr
\code{digits.I2} \tab 0 \tab number of digits for I-squared measure
  \cr
\code{digits.pval} \tab 3 \tab number of digits for p-values \cr
\code{CIbracket}, \tab "(" \tab \cr
\code{CIseparator} \tab "-" \tab print confidence intervals as
  "\code{(.-.)}" \cr
\code{zero.pval}, \tab FALSE \tab print p-values with leading zero
\cr
\code{JAMA.pval}, \tab TRUE \tab round p-values to three digits
  (for 0.001 < p \eqn{\le} 0.01) \cr
  \tab \tab or two digits (p > 0.01) \cr
\code{header.line}, \tab TRUE \tab print header line
}

IQWiG, General Methods 5 settings:
\tabular{lll}{
\bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
\code{method.random.ci} \tab "HK" \tab Hartung-Knapp method \cr
\code{prediction} \tab TRUE \tab Prediction interval \cr
}

IQWiG, General Methods 6 settings:
\tabular{lll}{
\bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
\code{method.random.ci} \tab "HK" \tab Hartung-Knapp method \cr
\code{adhoc.hakn.ci} \tab "IQWiG6" \tab \emph{ad hoc} variance correction \cr
\code{method.tau} \tab "PM" \tab Paule-Mandel estimator for
  between-study variance \cr
\code{prediction} \tab TRUE \tab Prediction interval \cr
}

Settings for gene expression data:
\tabular{lll}{
\bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
\code{scientific.pval} \tab TRUE \tab Scientific notation for p-values \cr
\code{method.tau.ci} \tab FALSE \tab
  no confidence interval for between-study \cr
 \tab \tab heterogeneity variance \cr
}

IVhet settings:
\tabular{lll}{
\bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
\code{method.common.ci} \tab "IVhet" \tab inverse variance heterogeneity \cr
\code{text.common} \tab "IVhet model" \tab  \cr
\code{text.w.common} \tab "IVhet" \tab  \cr
}

Settings for \bold{meta}, version 4 or below:
\tabular{lll}{
\bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
\code{method.tau} \tab "DL" \tab DerSimonian-Laird estimator \cr
\code{method.I2} \tab "Q" \tab Use Q to calculate I-squared \cr
\code{method.predict} \tab "HTS" \tab prediction interval with \emph{k-2}
  degrees \cr
  \tab \tab of freedom \cr
\code{exact.smd} \tab FALSE \tab Use exact formula for standardised mean \cr
  \tab \tab difference (White and Thomas, 2005) \cr
\code{text.common} \tab "Fixed effect model" \tab \cr
\code{text.w.common} \tab "fixed" \tab \cr
\code{warn.deprecated} \tab FALSE \tab Do not print warnings for deprecated
  \cr
  \tab \tab arguments
}

Settings for \bold{meta}, version 7.0-0 or below:
\tabular{lll}{
\bold{Argument} \tab \bold{Value} \tab \bold{Comment} \cr
\code{method.tau} \tab "REML" \tab REML estimator \cr
\code{method.I2} \tab "Q" \tab Use Q to calculate I-squared \cr
\code{method.predict} \tab "HTS" \tab prediction interval with \emph{k-2}
  degrees \cr
  \tab \tab of freedom \cr
\code{exact.smd} \tab TRUE \tab Use exact formula for standardised mean \cr
  \tab \tab difference (White and Thomas, 2005) \cr
\code{text.common} \tab "Common effect model" \tab \cr
\code{text.w.common} \tab "common" \tab \cr
\code{warn.deprecated} \tab FALSE \tab Do not print warnings for deprecated
  \cr
  \tab \tab arguments
}

A list of all arguments with current settings is printed using the
command \code{settings.meta()}.

In order to reset all settings of R package \bold{meta} the command
\code{settings.meta(reset = TRUE)} can be used.
}
\examples{
# Get listing of current settings
#
settings.meta()

# Meta-analyses using default settings
#
metabin(10, 20, 15, 20)
metaprop(4, 20)
metabin(10, 20, 15, 20, sm = "RD")
metaprop(4, 20, sm = "PLN")

# Change summary measure for R functions metabin and metaprop
# and store old settings
#
oldset <- settings.meta(smbin = "RD", smprop = "PLN")
#
metabin(10, 20, 15, 20)
metaprop(4, 20)

# Use old settings
#
settings.meta(oldset)

# Change level used to calculate confidence intervals
# (99\%-CI for studies, 99.9\%-CI for pooled effects)
#
metagen(1:3, 2:4 / 10, sm = "MD")
settings.meta(level = 0.99, level.ma = 0.999)
metagen(1:3, 2:4 / 10, sm = "MD")

# Always print a prediction interval
#
settings.meta(prediction = TRUE)
metagen(1:3, 2:4 / 10, sm = "MD")
metagen(4:6, 4:2 / 10, sm = "MD")

# Try to set unknown argument results in a warning
#
try(settings.meta(unknownarg = TRUE))

# Reset to default settings of R package meta
#
settings.meta("reset")
metabin(10, 20, 15, 20)
metaprop(4, 20)
metagen(1:3, 2:4 / 10, sm = "MD")

# Do not back transform results (e.g. print log odds ratios instead
# of odds ratios, print transformed correlations / proportions
# instead of correlations / proportions)
#
settings.meta(backtransf = FALSE)
metabin(10, 20, 15, 20)
metaprop(4, 20)
metacor(c(0.85, 0.7, 0.95), c(20, 40, 10))

# Forest plot using RevMan 5 style
#
settings.meta("RevMan5")
forest(metagen(1:3, 2:4 / 10, sm = "MD", common = FALSE),
  label.left = "Favours A", label.right = "Favours B",
  colgap.studlab = "2cm", colgap.forest.left = "0.2cm")

# Forest plot using JAMA style
#
settings.meta("JAMA")
forest(metagen(1:3, 2:4 / 10, sm = "MD", common = FALSE),
  label.left = "Favours A", label.right = "Favours B",
  colgap.studlab = "2cm", colgap.forest.left = "0.2cm")

# Use slightly different layout for confidence intervals
# (especially useful if upper confidence limit can be negative)
#
settings.meta(CIseparator = " - ")
forest(metagen(-(1:3), 2:4 / 10, sm = "MD", common = FALSE),
  label.left = "Favours A", label.right = "Favours B",
  colgap.studlab = "2cm", colgap.forest.left = "0.2cm")

# Use old settings
#
settings.meta(oldset)

}
\references{
Doi SAR, Barendregt JJ, Khan S, Thalib L, Williams GM (2015):
Advances in the meta-analysis of heterogeneous clinical trials I:
The inverse variance heterogeneity model.
\emph{Contemporary Clinical Trials},
\bold{45}, 130--8

White IR, Thomas J (2005):
Standardized mean differences in individually-randomized and
cluster-randomized trials, with applications to meta-analysis.
\emph{Clinical Trials},
\bold{2}, 141--51
}
\seealso{
\code{\link{gs}}, \code{\link{forest.meta}},
  \code{\link{print.meta}}, \code{\link{labels.meta}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
