% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaprop.R
\name{metaprop}
\alias{metaprop}
\title{Meta-analysis of single proportions}
\usage{
metaprop(
  event,
  n,
  studlab,
  data = NULL,
  subset = NULL,
  exclude = NULL,
  cluster = NULL,
  rho = 0,
  weights = NULL,
  weights.common = weights,
  weights.random = weights,
  method,
  sm = gs("smprop"),
  incr = gs("incr"),
  method.incr = gs("method.incr"),
  method.ci = gs("method.ci.prop"),
  level = gs("level"),
  common = gs("common"),
  random = gs("random") | !is.null(tau.preset),
  overall = common | random,
  overall.hetstat = if (is.null(gs("overall.hetstat"))) common | random else
    gs("overall.hetstat"),
  prediction = gs("prediction") | !missing(method.predict),
  method.tau = ifelse(!is.na(charmatch(tolower(method), "glmm", nomatch = NA)), "ML",
    gs("method.tau")),
  method.tau.ci = gs("method.tau.ci"),
  level.hetstat = gs("level.hetstat"),
  tau.preset = NULL,
  TE.tau = NULL,
  tau.common = gs("tau.common"),
  detail.tau = NULL,
  method.I2 = gs("method.I2"),
  level.ma = gs("level.ma"),
  method.common.ci = gs("method.common.ci"),
  method.random.ci = gs("method.random.ci"),
  adhoc.hakn.ci = gs("adhoc.hakn.ci"),
  level.predict = gs("level.predict"),
  method.predict = gs("method.predict"),
  adhoc.hakn.pi = gs("adhoc.hakn.pi"),
  seed.predict = NULL,
  null.effect = NA,
  method.bias = gs("method.bias"),
  backtransf = gs("backtransf"),
  pscale = 1,
  text.common = gs("text.common"),
  text.random = gs("text.random"),
  text.predict = gs("text.predict"),
  text.w.common = gs("text.w.common"),
  text.w.random = gs("text.w.random"),
  title = gs("title"),
  complab = gs("complab"),
  outclab = "",
  label.left = gs("label.left"),
  label.right = gs("label.right"),
  col.label.left = gs("col.label.left"),
  col.label.right = gs("col.label.right"),
  subgroup,
  subgroup.name = NULL,
  print.subgroup.name = gs("print.subgroup.name"),
  sep.subgroup = gs("sep.subgroup"),
  test.subgroup = gs("test.subgroup"),
  prediction.subgroup = gs("prediction.subgroup"),
  seed.predict.subgroup = NULL,
  byvar,
  hakn,
  adhoc.hakn,
  keepdata = gs("keepdata"),
  warn = gs("warn"),
  warn.deprecated = gs("warn.deprecated"),
  control = NULL,
  ...
)
}
\arguments{
\item{event}{Number of events.}

\item{n}{Number of observations.}

\item{studlab}{An optional vector with study labels.}

\item{data}{An optional data frame containing the study
information, i.e., event and n.}

\item{subset}{An optional vector specifying a subset of studies to
be used.}

\item{exclude}{An optional vector specifying studies to exclude
from meta-analysis, however, to include in printouts and forest
plots.}

\item{cluster}{An optional vector specifying which estimates come
from the same cluster resulting in the use of a three-level
meta-analysis model.}

\item{rho}{Assumed correlation of estimates within a cluster.}

\item{weights}{A single numeric or vector with user-specified weights.}

\item{weights.common}{User-specified weights (common effect model).}

\item{weights.random}{User-specified weights (random effects model).}

\item{method}{A character string indicating which method is to be
used for pooling of studies. One of \code{"Inverse"} and
\code{"GLMM"}, can be abbreviated.}

\item{sm}{A character string indicating which summary measure
(\code{"PLOGIT"}, \code{"PAS"}, \code{"PFT"}, \code{"PLN"}, or
\code{"PRAW"}) is to be used for pooling of studies, see Details.}

\item{incr}{A numeric which is added to event number and sample
size of studies with zero or all events, i.e., studies with an
event probability of either 0 or 1. Or a numeric vector with the
continuity correction for each study.}

\item{method.incr}{A character string indicating which continuity
correction method should be used (\code{"only0"},
\code{"if0all"}, or \code{"all"}), see Details.}

\item{method.ci}{A character string indicating which method is used
to calculate confidence intervals for individual studies, see
Details.}

\item{level}{The level used to calculate confidence intervals for
individual studies.}

\item{common}{A logical indicating whether a common effect
meta-analysis should be conducted.}

\item{random}{A logical indicating whether a random effects
meta-analysis should be conducted.}

\item{overall}{A logical indicating whether overall summaries
should be reported. This argument is useful in a meta-analysis
with subgroups if overall results should not be reported.}

\item{overall.hetstat}{A logical value indicating whether to print
heterogeneity measures for overall treatment comparisons. This
argument is useful in a meta-analysis with subgroups if
heterogeneity statistics should only be printed on subgroup
level.}

\item{prediction}{A logical indicating whether a prediction
interval should be printed.}

\item{method.tau}{A character string indicating which method is
used to estimate the between-study variance \eqn{\tau^2} and its
square root \eqn{\tau} (see \code{\link{meta-package}}).}

\item{method.tau.ci}{A character string indicating which method is
used to estimate the confidence interval of \eqn{\tau^2} and
\eqn{\tau} (see \code{\link{meta-package}}).}

\item{level.hetstat}{The level used to calculate confidence intervals
for heterogeneity statistics.}

\item{tau.preset}{Prespecified value for the square root of the
between-study variance \eqn{\tau^2}.}

\item{TE.tau}{Overall treatment effect used to estimate the
between-study variance tau-squared.}

\item{tau.common}{A logical indicating whether tau-squared should
be the same across subgroups.}

\item{detail.tau}{Detail on between-study variance estimate.}

\item{method.I2}{A character string indicating which method is
used to estimate the heterogeneity statistic I\eqn{^2}. Either
\code{"Q"} or \code{"tau2"}, can be abbreviated
(see \code{\link{meta-package}}).}

\item{level.ma}{The level used to calculate confidence intervals
for meta-analysis estimates.}

\item{method.common.ci}{A character string indicating which method
is used to calculate confidence interval and test statistic for
common effect estimate (see \code{\link{meta-package}}).}

\item{method.random.ci}{A character string indicating which method
is used to calculate confidence interval and test statistic for
random effects estimate (see \code{\link{meta-package}}).}

\item{adhoc.hakn.ci}{A character string indicating whether an
\emph{ad hoc} variance correction should be applied in the case
of an arbitrarily small Hartung-Knapp variance estimate (see
\code{\link{meta-package}}).}

\item{level.predict}{The level used to calculate prediction
interval for a new study.}

\item{method.predict}{A character string indicating which method is
used to calculate a prediction interval (see
\code{\link{meta-package}}).}

\item{adhoc.hakn.pi}{A character string indicating whether an
\emph{ad hoc} variance correction should be applied for
prediction interval (see \code{\link{meta-package}}).}

\item{seed.predict}{A numeric value used as seed to calculate
bootstrap prediction interval (see \code{\link{meta-package}}).}

\item{null.effect}{A numeric value specifying the effect under the
null hypothesis.}

\item{method.bias}{A character string indicating which test is to
be used. Either \code{"Begg"}, \code{"Egger"}, or
\code{"Thompson"}, can be abbreviated. See function
\code{\link{metabias}}.}

\item{backtransf}{A logical indicating whether results for
transformed proportions (argument \code{sm != "PRAW"}) should be
back transformed in printouts and plots. If TRUE (default),
results will be presented as proportions; otherwise transformed
proportions will be shown. See Details for presentation of
confidence intervals.}

\item{pscale}{A numeric defining a scaling factor for printing of
single event probabilities.}

\item{text.common}{A character string used in printouts and forest
plot to label the pooled common effect estimate.}

\item{text.random}{A character string used in printouts and forest
plot to label the pooled random effects estimate.}

\item{text.predict}{A character string used in printouts and forest
plot to label the prediction interval.}

\item{text.w.common}{A character string used to label weights of
common effect model.}

\item{text.w.random}{A character string used to label weights of
random effects model.}

\item{title}{Title of meta-analysis / systematic review.}

\item{complab}{Comparison label.}

\item{outclab}{Outcome label.}

\item{label.left}{Graph label on left side of null effect in forest plot.}

\item{label.right}{Graph label on right side of null effect in forest plot.}

\item{col.label.left}{The colour of the graph label on the left side of
the null effect.}

\item{col.label.right}{The colour of the graph label on the right side of
the null effect.}

\item{subgroup}{An optional vector to conduct a meta-analysis with
subgroups.}

\item{subgroup.name}{A character string with a name for the
subgroup variable.}

\item{print.subgroup.name}{A logical indicating whether the name of
the subgroup variable should be printed in front of the group
labels.}

\item{sep.subgroup}{A character string defining the separator
between name of subgroup variable and subgroup label.}

\item{test.subgroup}{A logical value indicating whether to print
results of test for subgroup differences.}

\item{prediction.subgroup}{A logical indicating whether prediction
intervals should be printed for subgroups.}

\item{seed.predict.subgroup}{A numeric vector providing seeds to
calculate bootstrap prediction intervals within subgroups. Must
be of same length as the number of subgroups.}

\item{byvar}{Deprecated argument (replaced by 'subgroup').}

\item{hakn}{Deprecated argument (replaced by 'method.random.ci').}

\item{adhoc.hakn}{Deprecated argument (replaced by 'adhoc.hakn.ci').}

\item{keepdata}{A logical indicating whether original data (set)
should be kept in meta object.}

\item{warn}{A logical indicating whether warnings should be printed
(e.g., if estimation problems exist in fitting a GLMM).}

\item{warn.deprecated}{A logical indicating whether warnings should
be printed if deprecated arguments are used.}

\item{control}{An optional list to control the iterative process to
estimate the between-study variance \eqn{\tau^2}. This argument
is passed on to \code{\link[metafor]{rma.uni}} or
\code{\link[metafor]{rma.glmm}}, respectively.}

\item{\dots}{Additional arguments passed on to
\code{\link[metafor]{rma.glmm}} function and to catch deprecated
arguments.}
}
\value{
An object of class \code{c("metaprop", "meta")} with corresponding
generic functions (see \code{\link{meta-object}}).
}
\description{
Calculation of an overall proportion from studies reporting a
single proportion. Inverse variance method and generalised linear
mixed model (GLMM) are available for pooling. For GLMMs, the
\code{\link[metafor]{rma.glmm}} function from R package
\bold{metafor} (Viechtbauer 2010) is called internally.
}
\details{
This function provides methods for common effect and random effects
meta-analysis of single proportions to calculate an overall
proportion. Note, you should use R function \code{\link{metabin}}
to compare proportions of pairwise comparisons instead of using
\code{metaprop} for each treatment arm separately which will break
randomisation in randomised controlled trials.

The following transformations of proportions are
implemented to calculate an overall proportion:

\itemize{
\item Logit transformation (\code{sm = "PLOGIT"}, default)
\item Arcsine transformation (\code{sm = "PAS"})
\item Freeman-Tukey Double arcsine transformation (\code{sm = "PFT"})
\item Log transformation (\code{sm = "PLN"})
\item No transformation (\code{sm = "PRAW"})
}

List elements \code{TE}, \code{TE.common}, \code{TE.random}, etc.,
contain the transformed proportions. In printouts and plots these
values are back transformed if argument \code{backtransf = TRUE}
(default).

A generalised linear mixed model (GLMM) - more specific, a random
intercept logistic regression model - can be utilised for the
meta-analysis of proportions (Stijnen et al., 2010). This is the
default method for the logit transformation (argument \code{sm =
"PLOGIT"}). Internally, the \code{\link[metafor]{rma.glmm}}
function from R package \bold{metafor} is called to fit a GLMM.

Classic meta-analysis (Borenstein et al., 2010) utilising the
(un)transformed proportions and corresponding standard errors in
the inverse variance method is conducted by calling the
\code{\link{metagen}} function internally. This is the only
available method for all transformations but the logit
transformation. The classic meta-analysis model with logit
transformed proportions is used by setting argument \code{method =
"Inverse"}.

A three-level random effects meta-analysis model (Van den Noortgate
et al., 2013) is utilised if argument \code{cluster} is used and at
least one cluster provides more than one estimate. Internally,
\code{\link[metafor]{rma.mv}} is called to conduct the analysis and
\code{\link[metafor]{weights.rma.mv}} with argument \code{type =
"rowsum"} is used to calculate random effects weights.

Default settings are utilised for several arguments (assignments
using \code{\link{gs}} function). These defaults can be changed for
the current R session using the \code{\link{settings.meta}}
function.

Furthermore, R function \code{\link{update.meta}} can be used to
rerun a meta-analysis with different settings.

\subsection{Choice of transformation / meta-analysis method}{

Contradictory recommendations on the use of transformations of
proportions have been published in the literature. For example,
Barendregt et al. (2013) recommend the use of the Freeman-Tukey
double arcsine transformation instead of the logit transformation
whereas Warton & Hui (2011) strongly advise to use generalised
linear mixed models with the logit transformation instead of the
arcsine transformation.

Schwarzer et al. (2019) describe seriously misleading results in a
meta-analysis with very different sample sizes due to problems with
the back-transformation of the Freeman-Tukey transformation which
requires a single sample size (Miller, 1978). Accordingly,
Schwarzer et al. (2019) also recommend to use GLMMs for the
meta-analysis of single proportions, however, admit that individual
study weights are not available with this method. Meta-analysts
which require individual study weights should consider the inverse
variance method with the arcsine or logit transformation.

In order to prevent misleading conclusions for the Freeman-Tukey
double arcsine transformation, sensitivity analyses using other
transformations or using a range of sample sizes should be
conducted (Schwarzer et al., 2019).
}

\subsection{Continuity correction}{

Three approaches are available to apply a continuity correction:
\itemize{
\item Only studies with a zero cell count (\code{method.incr =
  "only0"})
\item All studies if at least one study has a zero cell count
  (\code{method.incr = "if0all"})
\item All studies irrespective of zero cell counts
  (\code{method.incr = "all"})
}

If the summary measure is equal to "PLOGIT", "PLN", or "PRAW", the
continuity correction is applied if a study has either zero or all
events, i.e., an event probability of either 0 or 1.

By default, 0.5 is used as continuity correction (argument
\code{incr}). This continuity correction is used both to calculate
individual study results with confidence limits and to conduct
meta-analysis based on the inverse variance method. For GLMMs no
continuity correction is used. Furthermore, the value of \code{incr} is
only considered in the calculation of confidence intervals for individual
studies if \code{method.ci = "NAsm"} (see next subsection).
}

\subsection{Confidence intervals for individual studies}{

Various methods are available to calculate confidence intervals for
individual study results (see Agresti & Coull 1998 and Newcombe
1988):
\itemize{
\item Clopper-Pearson interval also called 'exact' binomial
  interval (\code{method.ci = "CP"}, default)
\item Wilson Score interval (\code{method.ci = "WS"})
\item Wilson Score interval with continuity correction
  (\code{method.ci = "WSCC"})
\item Agresti-Coull interval (\code{method.ci = "AC"})
\item Simple approximation interval (\code{method.ci = "SA"})
\item Simple approximation interval with continuity correction
  (\code{method.ci = "SACC"})
\item Normal approximation interval based on summary measure,
  i.e. defined by argument \code{sm} (\code{method.ci = "NAsm"})
}

Note, with exception of the normal approximation based on the
summary measure, i.e. \code{method.ci = "NAsm"}, the same
confidence interval is calculated for individual studies for any
summary measure (argument \code{sm}) as only number of events and
observations are used in the calculation disregarding the chosen
transformation. Furthermore, the continuity correction 

Results will be presented for transformed proportions if argument
\code{backtransf = FALSE}. In this case, argument \code{method.ci =
"NAsm"} is used, i.e. confidence intervals based on the normal
approximation based on the summary measure.
}

\subsection{Subgroup analysis}{

Argument \code{subgroup} can be used to conduct subgroup analysis for
a categorical covariate. The \code{\link{metareg}} function can be
used instead for more than one categorical covariate or continuous
covariates.
}

\subsection{Specify the null hypothesis of test for an overall proportion}{

Argument \code{null.effect} can be used to specify the proportion
used under the null hypothesis in a test for an overall effect.

By default (\code{null.effect = NA}), no hypothesis test is
conducted as it is unclear which value is a sensible choice for the
data at hand.  An overall proportion of 50\%, for example, could be
tested by setting argument \code{null.effect = 0.5}.

Note, all tests for an overall effect are two-sided with the
alternative hypothesis that the effect is unequal to
\code{null.effect}.
}

\subsection{Exclusion of studies from meta-analysis}{

Arguments \code{subset} and \code{exclude} can be used to exclude
studies from the meta-analysis. Studies are removed completely from
the meta-analysis using argument \code{subset}, while excluded
studies are shown in printouts and forest plots using argument
\code{exclude} (see Examples in \code{\link{metagen}}).
Meta-analysis results are the same for both arguments.
}

\subsection{Presentation of meta-analysis results}{

Internally, both common effect and random effects models are
calculated regardless of values choosen for arguments
\code{common} and \code{random}. Accordingly, the estimate
for the random effects model can be extracted from component
\code{TE.random} of an object of class \code{"meta"} even if
argument \code{random = FALSE}. However, all functions in R
package \bold{meta} will adequately consider the values for
\code{common} and \code{random}. E.g. function
\code{\link{print.meta}} will not print results for the random
effects model if \code{random = FALSE}.

Argument \code{pscale} can be used to rescale proportions, e.g.
\code{pscale = 1000} means that proportions are expressed as events
per 1000 observations. This is useful in situations with (very) low
event probabilities.

A prediction interval will only be shown if \code{prediction =
TRUE}.
}
}
\examples{
# Meta-analysis using generalised linear mixed model
#
metaprop(4:1, 10 * 1:4)

# Apply various classic meta-analysis methods to estimate
# proportions
#
m1 <- metaprop(4:1, 10 * 1:4, method = "Inverse")
m2 <- update(m1, sm = "PAS")
m3 <- update(m1, sm = "PRAW")
m4 <- update(m1, sm = "PLN")
m5 <- update(m1, sm = "PFT")
#
m1
m2
m3
m4
m5
#
forest(m1)
\dontrun{
forest(m2)
forest(m3)
forest(m3, pscale = 100)
forest(m4)
forest(m5)
}

# Do not back transform results, e.g. print logit transformed
# proportions if sm = "PLOGIT" and store old settings
#
oldset <- settings.meta(backtransf = FALSE)
#
m6  <- metaprop(4:1, c(10, 20, 30, 40), method = "Inverse")
m7  <- update(m6, sm = "PAS")
m8  <- update(m6, sm = "PRAW")
m9  <- update(m6, sm = "PLN")
m10 <- update(m6, sm = "PFT")
#
forest(m6)
\dontrun{
forest(m7)
forest(m8)
forest(m8, pscale = 100)
forest(m9)
forest(m10)
}

# Use old settings
#
settings.meta(oldset)

# Examples with zero events
#
m1 <- metaprop(c(0, 0, 10, 10), rep(100, 4), method = "Inverse")
m2 <- metaprop(c(0, 0, 10, 10), rep(100, 4), incr = 0.1, method = "Inverse")
#
m1
m2
#
\dontrun{
forest(m1)
forest(m2)
}

# Example from Miller (1978):
#
death <- c(3, 6, 10, 1)
animals <- c(11, 17, 21, 6)
#
m3 <- metaprop(death, animals, sm = "PFT")
forest(m3)

# Data examples from Newcombe (1998)
# - apply various methods to estimate confidence intervals for
#   individual studies
#
event <- c(81, 15, 0, 1)
n <- c(263, 148, 20, 29)
#
m1 <- metaprop(event, n, method.ci = "SA", method = "Inverse")
m2 <- update(m1, method.ci = "SACC")
m3 <- update(m1, method.ci = "WS")
m4 <- update(m1, method.ci = "WSCC")
m5 <- update(m1, method.ci = "CP")
#
lower <- round(logit2p(rbind(NA, m1$lower, m2$lower, NA, m3$lower,
  m4$lower, NA, m5$lower)), 4)
upper <- round(logit2p(rbind(NA, m1$upper, m2$upper, NA, m3$upper,
  m4$upper, NA, m5$upper)), 4)
#
tab1 <- data.frame(
  scen1 = meta:::formatCI(lower[, 1], upper[, 1]),
  scen2 = meta:::formatCI(lower[, 2], upper[, 2]),
  scen3 = meta:::formatCI(lower[, 3], upper[, 3]),
  scen4 = meta:::formatCI(lower[, 4], upper[, 4])
  )
names(tab1) <- c("r=81, n=263", "r=15, n=148",
  "r=0, n=20", "r=1, n=29")
row.names(tab1) <- c("Simple", "- SA", "- SACC",
  "Score", "- WS", "- WSCC", "Binomial", "- CP")
tab1[is.na(tab1)] <- ""
# Newcombe (1998), Table I, methods 1-5:
tab1

# Same confidence interval, i.e. unaffected by choice of summary
# measure
#
print(metaprop(event, n, method.ci = "WS", method = "Inverse"), ma = FALSE)
print(metaprop(event, n, sm = "PLN", method.ci = "WS"), ma = FALSE)
print(metaprop(event, n, sm = "PFT", method.ci = "WS"), ma = FALSE)
print(metaprop(event, n, sm = "PAS", method.ci = "WS"), ma = FALSE)
print(metaprop(event, n, sm = "PRAW", method.ci = "WS"), ma = FALSE)

# Different confidence intervals as argument sm = "NAsm"
#
print(metaprop(event, n, method.ci = "NAsm", method = "Inverse"), ma = FALSE)
print(metaprop(event, n, sm = "PLN", method.ci = "NAsm"), ma = FALSE)
print(metaprop(event, n, sm = "PFT", method.ci = "NAsm"), ma = FALSE)
print(metaprop(event, n, sm = "PAS", method.ci = "NAsm"), ma = FALSE)
print(metaprop(event, n, sm = "PRAW", method.ci = "NAsm"), ma = FALSE)

# Different confidence intervals as argument backtransf = FALSE.
# Accordingly, method.ci = "NAsm" used internally.
#
print(metaprop(event, n, method.ci = "WS", method = "Inverse"),
  ma = FALSE, backtransf = FALSE)
print(metaprop(event, n, sm = "PLN", method.ci = "WS"),
  ma = FALSE, backtransf = FALSE)
print(metaprop(event, n, sm = "PFT", method.ci = "WS"),
  ma = FALSE, backtransf = FALSE)
print(metaprop(event, n, sm = "PAS", method.ci = "WS"),
  ma = FALSE, backtransf = FALSE)
print(metaprop(event, n, sm = "PRAW", method.ci = "WS"),
  ma = FALSE, backtransf = FALSE)

# Same results (printed on original and log scale, respectively)
#
print(metaprop(event, n, sm = "PLN", method.ci = "NAsm"), ma = FALSE)
print(metaprop(event, n, sm = "PLN"), ma = FALSE, backtransf = FALSE)
# Results for first study (on log scale)
round(log(c(0.3079848, 0.2569522, 0.3691529)), 4)

# Print results as events per 1000 observations
#
print(metaprop(6:8, c(100, 1200, 1000), method = "Inverse"),
  pscale = 1000, digits = 1)

}
\references{
Agresti A & Coull BA (1998):
Approximate is better than "exact" for interval estimation of
binomial proportions.
\emph{The American Statistician},
\bold{52}, 119--26

Barendregt JJ, Doi SA, Lee YY, Norman RE, Vos T (2013):
Meta-analysis of prevalence.
\emph{Journal of Epidemiology and Community Health},
\bold{67}, 974--8

Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2010):
A basic introduction to fixed-effect and random-effects models for
meta-analysis.
\emph{Research Synthesis Methods},
\bold{1}, 97--111

Freeman MF & Tukey JW (1950):
Transformations related to the angular and the square root.
\emph{Annals of Mathematical Statistics},
\bold{21}, 607--11

Miller JJ (1978):
The inverse of the Freeman-Tukey double arcsine transformation.
\emph{The American Statistician},
\bold{32}, 138

Newcombe RG (1998):
Two-sided confidence intervals for the single proportion:
comparison of seven methods.
\emph{Statistics in Medicine},
\bold{17}, 857--72

Pettigrew HM, Gart JJ, Thomas DG (1986):
The bias and higher cumulants of the logarithm of a binomial
variate.
\emph{Biometrika},
\bold{73}, 425--35

Schwarzer G, Chemaitelly H, Abu-Raddad LJ, Rücker G (2019):
Seriously misleading results using inverse of Freeman-Tukey double
arcsine transformation in meta-analysis of single proportions.
\emph{Research Synthesis Methods},
\bold{10}, 476--83

Stijnen T, Hamza TH, Ozdemir P (2010):
Random effects meta-analysis of event outcome in the framework of
the generalized linear mixed model with applications in sparse
data.
\emph{Statistics in Medicine},
\bold{29}, 3046--67

Van den Noortgate W, López-López JA, Marín-Martínez F, Sánchez-Meca J (2013):
Three-level meta-analysis of dependent effect sizes.
\emph{Behavior Research Methods},
\bold{45}, 576--94

Viechtbauer W (2010):
Conducting meta-analyses in R with the metafor package.
\emph{Journal of Statistical Software},
\bold{36}, 1--48

Warton DI, Hui FKC (2011):
The arcsine is asinine: the analysis of proportions in ecology.
\emph{Ecology},
\bold{92}, 3--10
}
\seealso{
\code{\link{meta-package}}, \code{\link{update.meta}},
  \code{\link{metacont}}, \code{\link{metagen}},
  \code{\link{print.meta}}, \code{\link{forest.meta}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
