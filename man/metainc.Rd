% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metainc.R
\name{metainc}
\alias{metainc}
\title{Meta-analysis of incidence rates}
\usage{
metainc(
  event.e,
  time.e,
  event.c,
  time.c,
  studlab,
  data = NULL,
  subset = NULL,
  exclude = NULL,
  cluster = NULL,
  rho = 0,
  weights = NULL,
  weights.common = weights,
  weights.random = weights,
  method = if (sm == "IRSD") "Inverse" else "MH",
  sm = gs("sminc"),
  incr = gs("incr"),
  method.incr = gs("method.incr"),
  incr.e = if (length(incr) > 1) incr else NULL,
  incr.c = if (length(incr) > 1) incr else NULL,
  model.glmm = "UM.FS",
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
  method.bias = gs("method.bias"),
  n.e = NULL,
  n.c = NULL,
  backtransf = if (sm == "IRSD") FALSE else gs("backtransf"),
  irscale = 1,
  irunit = "person-years",
  text.common = gs("text.common"),
  text.random = gs("text.random"),
  text.predict = gs("text.predict"),
  text.w.common = gs("text.w.common"),
  text.w.random = gs("text.w.random"),
  title = gs("title"),
  complab = gs("complab"),
  outclab = "",
  label.e = gs("label.e"),
  label.c = gs("label.c"),
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
\item{event.e}{Number of events in experimental group or an R object
created with \code{\link{pairwise}}.}

\item{time.e}{Person time at risk in experimental group.}

\item{event.c}{Number of events in control group.}

\item{time.c}{Person time at risk in control group.}

\item{studlab}{An optional vector with study labels.}

\item{data}{An optional data frame containing the study
information, i.e., event.e, time.e, event.c, and time.c.}

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
used for pooling of studies. One of \code{"MH"},
\code{"Inverse"}, \code{"Cochran"}, or \code{"GLMM"} can be
abbreviated.}

\item{sm}{A character string indicating which summary measure
(\code{"IRR"}, \code{"IRD"}, \code{"IRSD"}, or \code{"VE"}) is to
be used for pooling of studies, see Details.}

\item{incr}{A numerical value which is added to cell frequencies
for studies with a zero cell count or a numeric vector with the continuity
correction for each study, see Details.}

\item{method.incr}{A character string indicating which continuity
correction method should be used (\code{"only0"},
\code{"if0all"}, \code{"all"}, or \code{"user"}), see Details.}

\item{incr.e}{Continuity correction in experimental group, see Details.}

\item{incr.c}{Continuity correction in control group, see Details.}

\item{model.glmm}{A character string indicating which GLMM should
be used. One of \code{"UM.FS"}, \code{"UM.RS"}, and
\code{"CM.EL"}, see Details.}

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

\item{method.bias}{A character string indicating which test is to
be used. Either \code{"Begg"}, \code{"Egger"}, or
\code{"Thompson"}, can be abbreviated. See function
\code{\link{metabias}}.}

\item{n.e}{Number of observations in experimental group (optional).}

\item{n.c}{Number of observations in control group (optional).}

\item{backtransf}{A logical indicating whether results for
incidence rate ratio (\code{sm = "IRR"}) and vaccine efficacy or
vaccine effectiveness (\code{sm = "VE"}) should be back
transformed in printouts and plots. If TRUE (default), results
will be presented as incidence rate ratios or vaccine efficacy /
effectiveness; otherwise log incidence rate ratios or log vaccine
rate ratios will be shown.}

\item{irscale}{A numeric defining a scaling factor for printing of
incidence rate differences.}

\item{irunit}{A character string specifying the time unit used to
calculate rates, e.g. person-years.}

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

\item{label.e}{Label for experimental group.}

\item{label.c}{Label for control group.}

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

\item{adhoc.hakn}{Deprecated argument (replaced by
'adhoc.hakn.ci').}

\item{keepdata}{A logical indicating whether original data (set)
should be kept in meta object.}

\item{warn}{A logical indicating whether warnings should be printed
(e.g., if \code{incr} is added to studies with zero cell
frequencies or if estimation problems exist in fitting a GLMM).}

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
An object of class \code{c("metainc", "meta")} with corresponding
generic functions (see \code{\link{meta-object}}).
}
\description{
Calculation of common effect and random effects estimates (incidence
rate ratio or incidence rate difference) for meta-analyses with
event counts.  Mantel-Haenszel, Cochran, inverse variance method,
and generalised linear mixed model (GLMM) are available for
pooling. For GLMMs, the \code{\link[metafor]{rma.glmm}} function
from R package \bold{metafor} (Viechtbauer 2010) is called
internally.
}
\details{
Calculation of common and random effects estimates for meta-analyses
comparing two incidence rates.

The following measures of treatment effect are available:
\itemize{
\item Incidence Rate Ratio (\code{sm = "IRR"})
\item Incidence Rate Difference (\code{sm = "IRD"})
\item Square root transformed Incidence Rate Difference (\code{sm =
  "IRSD"})
\item Vaccine efficacy or vaccine effectiveness (\code{sm = "VE"})
}

Note, log incidence rate ratio (logIRR) and log vaccine ratio
(logVR) are mathematical identical, however, back-transformed
results differ as vaccine efficacy or effectiveness is defined as
\code{VE = 100 * (1 - IRR)}.

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

\subsection{Meta-analysis method}{

By default, both common effect and random effects models are
considered (see arguments \code{common} and \code{random}). If
\code{method} is \code{"MH"} (default), the Mantel-Haenszel method
is used to calculate the common effect estimate (Greenland &
Robbins, 1985); if \code{method} is \code{"Inverse"}, inverse
variance weighting is used for pooling; if \code{method} is
\code{"Cochran"}, the Cochran method is used for pooling
(Bayne-Jones, 1964, Chapter 8). For these three methods, the random
effects estimate is always based on the inverse variance method.

A distinctive and frequently overlooked advantage of incidence
rates is that individual patient data (IPD) can be extracted from
count data. Accordingly, statistical methods for IPD, i.e.,
generalised linear mixed models, can be utilised in a meta-analysis
of incidence rate ratios (Stijnen et al., 2010). These methods are
available (argument \code{method = "GLMM"}) for the common effect
and random effects model by calling the
\code{\link[metafor]{rma.glmm}} function from R package
\bold{metafor} internally.

Three different GLMMs are available for meta-analysis of incidence
rate ratios using argument \code{model.glmm} (which corresponds to
argument \code{model} in the \code{\link[metafor]{rma.glmm}}
function):
\tabular{cl}{
1. \tab Poisson regression model with fixed study effects (default)
 \cr
 \tab (\code{model.glmm = "UM.FS"}, i.e., \bold{U}nconditional
 \bold{M}odel - \bold{F}ixed \bold{S}tudy effects) \cr
2. \tab Mixed-effects Poisson regression model with random study
 effects \cr
 \tab (\code{model.glmm = "UM.RS"}, i.e., \bold{U}nconditional
 \bold{M}odel - \bold{R}andom \bold{S}tudy effects) \cr
3. \tab Generalised linear mixed model (conditional Poisson-Normal)
 \cr
 \tab (\code{model.glmm = "CM.EL"}, i.e., \bold{C}onditional
  \bold{M}odel - \bold{E}xact \bold{L}ikelihood)
}

Details on these three GLMMs as well as additional arguments which
can be provided using argument '\code{\dots}' in \code{metainc}
are described in \code{\link[metafor]{rma.glmm}} where you can also
find information on the iterative algorithms used for estimation.
Note, regardless of which value is used for argument
\code{model.glmm}, results for two different GLMMs are calculated:
common effect model (with fixed treatment effect) and random effects
model (with random treatment effects).
}

\subsection{Continuity correction}{

Four approaches are available to apply a continuity correction:
\itemize{
\item Only studies with a zero cell count (\code{method.incr =
  "only0", default})
\item All studies if at least one study has a zero cell count
  (\code{method.incr = "if0all"})
\item All studies irrespective of zero cell counts
  (\code{method.incr = "all"})
\item Use values provided in arguments \code{incr.e} and \code{incr.c}
  (\code{method.incr = "user"})
}

For studies with a zero cell count, by default, 0.5 is added to all
cell frequencies of these studies (argument \code{incr}). This
continuity correction is used both to calculate individual study
results with confidence limits and to conduct meta-analysis based
on the inverse variance method. For Mantel-Haenszel method, Cochran
method, and GLMMs, nothing is added to zero cell counts.
Accordingly, estimates for these methods are not defined if the
number of events is zero in all studies either in the experimental
or control group.
}

\subsection{Subgroup analysis}{

Argument \code{subgroup} can be used to conduct subgroup analysis for
a categorical covariate. The \code{\link{metareg}} function can be
used instead for more than one categorical covariate or continuous
covariates.
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

A prediction interval will only be shown if \code{prediction =
TRUE}.
}
}
\examples{
data(smoking)
m1 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
  data = smoking, studlab = study)
print(m1, digits = 2)

m2 <- update(m1, method = "Cochran")
print(m2, digits = 2)

data(lungcancer)
m3 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
  data = lungcancer, studlab = study)
print(m3, digits = 2)

# Redo Cochran meta-analysis with inflated standard errors
#
# All cause mortality
#
TEa <- log((smoking$d.smokers/smoking$py.smokers) /
  (smoking$d.nonsmokers/smoking$py.nonsmokers))
seTEa <- sqrt(1 / smoking$d.smokers + 1 / smoking$d.nonsmokers +
  2.5 / smoking$d.nonsmokers)
metagen(TEa, seTEa, sm = "IRR", studlab = smoking$study)

# Lung cancer mortality
#
TEl <- log((lungcancer$d.smokers/lungcancer$py.smokers) /
  (lungcancer$d.nonsmokers/lungcancer$py.nonsmokers))
seTEl <- sqrt(1 / lungcancer$d.smokers + 1 / lungcancer$d.nonsmokers +
  2.25 / lungcancer$d.nonsmokers)
metagen(TEl, seTEl, sm = "IRR", studlab = lungcancer$study)

\dontrun{
# Meta-analysis using generalised linear mixed models
# (only if R packages 'metafor' and 'lme4' are available)

# Poisson regression model (fixed study effects)
#
m4 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
  data = smoking, studlab = study, method = "GLMM")
m4

# Mixed-effects Poisson regression model (random study effects)
#
update(m4, model.glmm = "UM.RS", nAGQ = 1)
#
# Generalised linear mixed model (conditional Poisson-Normal)
#
update(m4, model.glmm = "CM.EL")
}

}
\references{
Bayne-Jones S et al. (1964):
Smoking and Health: Report of the Advisory Committee to the Surgeon
General of the United States.
U-23 Department of Health, Education, and Welfare.
Public Health Service Publication No. 1103.

Greenland S & Robins JM (1985):
Estimation of a common effect parameter from sparse follow-up data.
\emph{Biometrics},
\bold{41}, 55--68

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
Conducting Meta-Analyses in R with the Metafor Package.
\emph{Journal of Statistical Software},
\bold{36}, 1--48
}
\seealso{
\code{\link{meta-package}}, \code{\link{metabin}},
  \code{\link{update.meta}}, \code{\link{print.meta}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
