% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta-package.R
\docType{package}
\name{meta-package}
\alias{meta-package}
\alias{meta}
\title{meta: Brief overview of methods and general hints}
\description{
R package \bold{meta} is a user-friendly general package providing
standard methods for meta-analysis and supporting Schwarzer et
al. (2015),
\url{https://link.springer.com/book/10.1007/978-3-319-21416-0}.
}
\details{
R package \bold{meta} (Schwarzer, 2007; Balduzzi et al., 2019)
provides the following statistical methods for meta-analysis.
\enumerate{
\item Common effect (also called fixed effect) and random effects model:
\itemize{
 \item Meta-analysis of continuous outcome data (\code{\link{metacont}})
 \item Meta-analysis of binary outcome data (\code{\link{metabin}})
 \item Meta-analysis of incidence rates (\code{\link{metainc}})
 \item Generic inverse variance meta-analysis (\code{\link{metagen}})
 \item Meta-analysis of single correlations (\code{\link{metacor}})
 \item Meta-analysis of single means (\code{\link{metamean}})
 \item Meta-analysis of single proportions (\code{\link{metaprop}})
 \item Meta-analysis of single incidence rates (\code{\link{metarate}})
}
\item Several plots for meta-analysis:
\itemize{
 \item Forest plot (\code{\link{forest.meta}}, \code{\link{forest.metabind}})
 \item Funnel plot (\code{\link{funnel.meta}})
 \item Galbraith plot / radial plot (\code{\link{radial.meta}})
 \item L'Abbe plot for meta-analysis with binary outcome data
  (\code{\link{labbe.metabin}}, \code{\link{labbe.default}})
 \item Baujat plot to explore heterogeneity in meta-analysis
  (\code{\link{baujat.meta}})
 \item Bubble plot to display the result of a meta-regression
  (\code{\link{bubble.metareg}})
}
\item Three-level meta-analysis model (Van den Noortgate et al.,
  2013)
\item Generalised linear mixed models (GLMMs) for binary and count
  data (Stijnen et al., 2010) (\code{\link{metabin}},
  \code{\link{metainc}}, \code{\link{metaprop}}, and
  \code{\link{metarate}})
\item Logistic regression with penalised likelihood for rare binary data
  (Evrenoglou et al., 2022) (\code{\link{metabin}},
\item Various estimators for the between-study variance
 \eqn{\tau^2} in a random effects model (Veroniki et al., 2016);
 see description of argument \code{method.tau} below
\item Two methods to estimate the I-squared statistic
  (Higgins and Thompson, 2002); see description of argument
  \code{method.I2} below
\item Hartung-Knapp method for random effects meta-analysis
 (Hartung & Knapp, 2001a,b), see description of arguments
 \code{method.random.ci} and \code{adhoc.hakn.ci} below
\item Kenward-Roger method for random effects meta-analysis
 (Partlett and Riley, 2017), see description of arguments
 \code{method.random.ci} and \code{method.predict} below
\item Inverse variance heterogeneity method (Doi et al., 2015),
 see description of argument \code{method.common.ci} below
\item Prediction interval for the treatment effect of a new study
 (Veroniki et al., 2019; Higgins et al., 2009; Partlett and Riley, 2017;
 Nagashima et al., 2019), see description of argument \code{method.predict}
 below
\item Statistical tests for funnel plot asymmetry
 (\code{\link{metabias.meta}}, \code{\link{metabias.rm5}}) and
 trim-and-fill method (\code{\link{trimfill.meta}},
 \code{\link{trimfill.default}}) to evaluate bias in meta-analysis
\item Meta-regression (\code{\link{metareg}})
\item Cumulative meta-analysis (\code{\link{metacum}}) and
  leave-one-out meta-analysis (\code{\link{metainf}})
\item Import data from RevMan Web (\code{\link{read.cdir}}), RevMan
  5 (\code{\link{read.rm5}}), see also \code{\link{metacr}} to
  conduct meta-analysis for a single comparison and outcome from a
  Cochrane review
}

R package \bold{meta} provides two vignettes:
\itemize{
 \item \code{vignette("meta-workflow")} with an overview of main
   functions,
 \item \code{vignette("meta-tutorial")} with up-to-date commands for
   Balduzzi et al. (2019).
}

Additional statistical meta-analysis methods are provided by add-on
R packages:
\itemize{
\item Frequentist methods for network meta-analysis (R package
  \bold{netmeta})
\item Statistical methods for sensitivity analysis in meta-analysis
  (R package \bold{metasens})
\item Statistical methods for meta-analysis of diagnostic accuracy
  studies with several cutpoints (R package \bold{diagmeta})
}

In the following, more details on available and default statistical
meta-analysis methods are provided and R function
\code{\link{settings.meta}} is briefly described which can be used
to change the default settings. Additional information on
meta-analysis objects and available summary measures can be found
on the help pages \code{\link{meta-object}} and
\code{\link{meta-sm}}.

\subsection{Estimation of between-study variance}{

The following methods are available in all meta-analysis functions
to estimate the between-study variance \eqn{\tau^2}.
\tabular{ll}{
\bold{Argument} \tab \bold{Method} \cr 
\code{method.tau = "REML"}
 \tab Restricted maximum-likelihood estimator (Viechtbauer, 2005) \cr
\tab (default) \cr
\code{method.tau = "PM"}
 \tab Paule-Mandel estimator (Paule and Mandel, 1982) \cr
\code{method.tau = "DL"}
 \tab DerSimonian-Laird estimator (DerSimonian and Laird, 1986) \cr
\code{method.tau = "ML"}
 \tab Maximum-likelihood estimator (Viechtbauer, 2005) \cr
\code{method.tau = "HS"}
 \tab Hunter-Schmidt estimator (Hunter and Schmidt, 2015) \cr
\code{method.tau = "SJ"}
 \tab Sidik-Jonkman estimator (Sidik and Jonkman, 2005) \cr
\code{method.tau = "HE"}
 \tab Hedges estimator (Hedges and Olkin, 1985) \cr
\code{method.tau = "EB"}
 \tab Empirical Bayes estimator (Morris, 1983)
}

For GLMMs, only the maximum-likelihood method is available.

Historically, the DerSimonian-Laird method was the de facto
standard to estimate the between-study variance \eqn{\tau^2} and is
the default in some software packages including Review Manager 5
(RevMan 5) and R package \bold{meta}, version 4 and below. However,
its role has been challenged and especially the REML and
Paule-Mandel estimators have been recommended (Veroniki et al.,
2016; Langan et al., 2019). Accordingly, the currenct default in R
package \bold{meta} is the REML estimator.

The following R command could be used to employ the Paule-Mandel
instead of the REML estimator in all meta-analyses of the current R
session:
\itemize{
\item \code{settings.meta(method.tau = "PM")}
}

Other estimators for \eqn{\tau^2} could be selected in a similar
way.

Note, for binary outcomes, two variants of the DerSimonian-Laird
estimator are available if the Mantel-Haenszel method is used for
pooling. If argument \code{Q.Cochrane = TRUE} (default), the
heterogeneity statistic Q is based on the Mantel-Haenszel instead
of the inverse variance estimator under the common effect
model. This is the estimator for \eqn{\tau^2} implemented in RevMan
5.
}


\subsection{Estimation of I-squared statistic}{

The following methods are available in all meta-analysis functions
to estimate the I-squared statistic (Higgins and Thompson, 2002).
\tabular{ll}{
\bold{Argument} \tab \bold{Method} \cr 
\code{method.I2 = "Q"}
 \tab Based on heterogeneity statistic Q (default) \cr
\code{method.I2 = "tau2"}
 \tab Based on between-study variance \eqn{\tau^2}
}

Using \code{method.I2 = "Q"} (Higgins and Thompson, 2002, section 3.3), the
value of I\eqn{^2} does not change if the estimate of \eqn{\tau^2} changes.
Furthermore, the value of I\eqn{^2} and the test of heterogeneity based on
the Q statistic are in agreement. R package \bold{metafor} uses the second
method (\code{method.I2 = "tau2"}) which is described in Higgins and Thompson
(2002), section 3.2. This method is more general in the way that the value
of I\eqn{^2} changes with the estimate of \eqn{\tau^2}.
}

\subsection{Confidence interval for common effect estimate}{

The following methods are available in all meta-analysis functions
to calculate a confidence interval for the common effect estimate.
\tabular{ll}{
\bold{Argument} \tab \bold{Method} \cr
\code{method.common.ci = "classic"} \tab Based on standard normal
  quantile \cr
\code{method.common.ci = "IVhet"} \tab Method by Doi et al. (2015) \cr
}

The inverse variance heterogeneity method by Doi et al. (2015) is only
available if argument \code{method = "Inverse"}.
}

\subsection{Confidence interval for random effects estimate}{

The following methods are available in all meta-analysis functions
to calculate a confidence interval for the random effects estimate.
\tabular{ll}{
\bold{Argument} \tab \bold{Method} \cr
\code{method.random.ci = "classic"} \tab Based on standard normal
  quantile \cr
\tab (DerSimonian and Laird, 1986) (default) \cr
\code{method.random.ci = "HK"} \tab Method by Hartung and Knapp
  (2001a/b) \cr
\code{method.random.ci = "KR"} \tab Kenward-Roger method (Partlett and
  Riley, 2017)
}

DerSimonian and Laird (1986) introduced the classic random effects
model using a quantile of the standard normal distribution to
calculate a confidence interval for the random effects
estimate. This method implicitly assumes that the weights in the
random effects meta-analysis are not estimated but
given. Particularly, the uncertainty in the estimation of the
between-study variance \eqn{\tau^2} is ignored.

Hartung and Knapp (2001a,b) proposed an alternative method for
random effects meta-analysis based on a refined variance estimator
for the treatment estimate and a quantile of a
\emph{t}-distribution with \emph{k-1} degrees of freedom where
\emph{k} corresponds to the number of studies in the
meta-analysis.

The Kenward-Roger method is only available for the REML estimator
(\code{method.tau = "REML"}) of the between-study variance
\eqn{\tau^2} (Partlett and Riley, 2017). This method is based on an
adjusted variance estimate for the random effects
estimate. Furthermore, a quantile of a \emph{t}-distribution with
adequately modified degrees of freedom is used to calculate the
confidence interval.

For GLMMs and three-level models, the Kenward-Roger method is not
available, but a method similar to Knapp and Hartung (2003) is used
if \code{method.random.ci = "HK"}. For this method, the variance
estimator is not modified, however, a quantile of a
\emph{t}-distribution with \emph{k-1} degrees of freedom is used;
see description of argument \code{test} in
\code{\link[metafor]{rma.glmm}} and \code{\link[metafor]{rma.mv}}.

Simulation studies (Hartung and Knapp, 2001a,b; IntHout et al.,
2014; Langan et al., 2019) show improved coverage probabilities of
the Hartung-Knapp method compared to the classic random effects
method. However, in rare settings with very homogeneous treatment
estimates, the Hartung-Knapp variance estimate can be arbitrarily
small resulting in a very narrow confidence interval (Knapp and
Hartung, 2003; Wiksten et al., 2016). In such cases, an \emph{ad
hoc} variance correction has been proposed by utilising the
variance estimate from the classic random effects model with the
Hartung-Knapp method (Knapp and Hartung, 2003; IQWiQ, 2022). An
alternative \emph{ad hoc} approach is to use the confidence
interval of the classic common or random effects meta-analysis if
it is wider than the interval from the Hartung-Knapp method
(Wiksten et al., 2016; Jackson et al., 2017).

Argument \code{adhoc.hakn.ci} can be used to choose the \emph{ad
hoc} correction for the Hartung-Knapp (HK) method:
\tabular{ll}{
\bold{Argument} \tab \bold{\emph{Ad hoc} method} \cr
\code{adhoc.hakn.ci = ""} \tab no \emph{ad hoc} correction (default)
  \cr
\code{adhoc.hakn.ci = "se"} \tab use variance correction if HK standard
 error is smaller \cr
 \tab than standard error from classic random effects
 \cr
 \tab meta-analysis (Knapp and Hartung, 2003) \cr
\code{adhoc.hakn.ci = "IQWiG6"} \tab use variance correction if HK
 confidence interval \cr
 \tab is narrower than CI from classic random effects model \cr
 \tab with DerSimonian-Laird estimator (IQWiG, 2022) \cr
\code{adhoc.hakn.ci = "ci"} \tab use wider confidence interval of
 classic random effects \cr
 \tab and HK meta-analysis \cr
 \tab (Hybrid method 2 in Jackson et al., 2017)
}

For GLMMs and three-level models, the \emph{ad hoc} variance
corrections are not available.
}

\subsection{Prediction interval}{

The following methods are available in all meta-analysis functions
to calculate a prediction interval for the treatment effect in a
single new study.
\tabular{ll}{
\bold{Argument} \tab \bold{Method} \cr
\code{method.predict = "V"} \tab Based on \emph{t}-distribution
  with \emph{k-1} degrees of freedom \cr
  \tab (Veroniki et al., 2019) (default) \cr
\code{method.predict = "HTS"} \tab Based on \emph{t}-distribution
  with \emph{k-2} degrees of freedom \cr
  \tab (Higgins et al., 2009) \cr
\code{method.predict = "HK"} \tab Based on Hartung-Knapp standard error
  and \cr
  \tab \emph{t}-distribution with \emph{k-1} degrees of freedom \cr
\code{method.predict = "HK-PR"} \tab Based on Hartung-Knapp standard error
  and \cr
  \tab \emph{t}-distribution with \emph{k-2} degrees of freedom \cr
  \tab (Partlett and Riley, 2017) \cr
\code{method.predict = "KR"} \tab Based on Kenward-Roger standard error
  and \cr
  \tab \emph{t}-distribution with approximate Kenward-Roger \cr
  \tab degrees of freedom \cr
\code{method.predict = "KR-PR"} \tab Based on Kenward-Roger standard error
  and \cr
  \tab \emph{t}-distribution with approximate Kenward-Roger \cr
  \tab degrees of freedom minus 1 (Partlett and Riley, 2017) \cr
\code{method.predict = "NNF"} \tab Bootstrap approach (Nagashima et
  al., 2019) \cr
\code{method.predict = "S"} \tab Based on standard normal quantile
  (Skipka, 2006)
}

By default (\code{method.predict = "V"}), the prediction interval
is based on a \emph{t}-distribution with \emph{k-1} degrees of
freedom where \emph{k} corresponds to the number of studies in the
meta-analysis (Veroniki et al., 2019). The method by Higgins et al., (2009),
which is based on a \emph{t}-distribution with \emph{k-2} degrees of freedom,
has been the default in R package \bold{meta}, version 7.0-0 or lower.

The Hartung-Knapp prediction intervals are also based on a
\emph{t}-distribution, however, use a different standard error.

The Kenward-Roger method is only available for the REML estimator
(\code{method.tau = "REML"}) of the between-study variance
\eqn{\tau^2} (Partlett and Riley, 2017). This method is based on an
adjusted variance estimate for the random effects
estimate. Furthermore, a quantile of a \emph{t}-distribution with
adequately modified degrees of freedom is used to calculate the prediction
interval.

The bootstrap approach is only available if R package \bold{pimeta}
is installed (Nagashima et al., 2019). Internally, the
\code{\link[pimeta]{pima}} function is called with argument
\code{method = "boot"}. Argument \code{seed.predict} can be used to
get a reproducible bootstrap prediction interval and argument
\code{seed.predict.subgroup} for reproducible bootstrap prediction
intervals in subgroups.

The method of Skipka (2006) ignores the uncertainty in the
estimation of the between-study variance \eqn{\tau^2} and thus has
too narrow limits for meta-analyses with a small number of studies.

For GLMMs and three-level models, only the methods by Veroniki et al.
(2019), Higgins et al. (2009) and Skipka (2006) are available. Argument
\code{method.predict = "V"} in R package \bold{meta} gives the same
prediction intervals as R functions \code{\link[metafor]{rma.glmm}} or
\code{\link[metafor]{rma.mv}} with argument \code{test = "t"}.

Note, in R package \bold{meta}, version 7.0-0 or lower, the methods
\code{method.predict = "HK-PR"} and \code{method.predict = "KR-PR"} have
been available as \code{method.predict = "HK"} and
\code{method.predict = "KR"}.

Argument \code{adhoc.hakn.pi} can be used to choose the \emph{ad
hoc} correction for the Hartung-Knapp method:

\tabular{ll}{
\bold{Argument} \tab \bold{\emph{Ad hoc} method} \cr
\code{adhoc.hakn.pi = ""} \tab no \emph{ad hoc} correction (default)
  \cr
\code{adhoc.hakn.pi = "se"} \tab use variance correction if HK
 standard error is smaller
}
}

\subsection{Confidence interval for the between-study variance}{

The following methods are available in all meta-analysis functions
to calculate a confidence interval for \eqn{\tau^2} and \eqn{\tau}.
\tabular{ll}{
\bold{Argument} \tab \bold{Method} \cr 
\code{method.tau.ci = "J"} \tab Method by Jackson (2013) \cr
\code{method.tau.ci = "BJ"} \tab Method by Biggerstaff and Jackson (2008) \cr
\code{method.tau.ci = "QP"} \tab Q-Profile method (Viechtbauer, 2007) \cr
\code{method.tau.ci = "PL"} \tab Profile-Likelihood method for three-level \cr
 \tab meta-analysis model (Van den Noortgate et al., 2013) \cr
\code{method.tau.ci = ""} \tab No confidence interval
}
The first three methods have been recommended by Veroniki et
al. (2016). By default, the Jackson method is used for the
DerSimonian-Laird estimator of \eqn{\tau^2} and the Q-profile
method for all other estimators of \eqn{\tau^2}.

The Profile-Likelihood method is the only method available for the
three-level meta-analysis model.

For GLMMs, no confidence intervals for \eqn{\tau^2} and \eqn{\tau}
are calculated.
}

\subsection{Change default settings for R session}{

R function \code{\link{settings.meta}} can be used to change the
previously described and several other default settings for the
current R session.

Some pre-defined general settings are available:
\itemize{
\item \code{settings.meta("RevMan5")}
\item \code{settings.meta("JAMA")}
\item \code{settings.meta("BMJ")}
\item \code{settings.meta("IQWiG5")}
\item \code{settings.meta("IQWiG6")}
\item \code{settings.meta("geneexpr")}
}

The first command can be used to reproduce meta-analyses from
Cochrane reviews conducted with \emph{Review Manager 5} (RevMan 5)
and specifies to use a RevMan 5 layout in forest plots.

The second command can be used to generate forest plots following
instructions for authors of the \emph{Journal of the American
Medical Association}. Study labels according to JAMA guidelines can be
generated using \code{\link{labels.meta}}.

The third command can be used to generate forest plots in the current layout
of the \emph{British Medical Journal}.

The next two commands implement the recommendations of the
Institute for Quality and Efficiency in Health Care (IQWiG),
Germany accordinging to General Methods 5 and 6, respectively
(\url{https://www.iqwig.de/en/about-us/methods/methods-paper/}).

The last setting can be used to print p-values in scientific
notation and to suppress the calculation of confidence intervals
for the between-study variance.

See \code{\link{settings.meta}} for more details on these
pre-defined general settings.

In addition, \code{\link{settings.meta}} can be used to define
individual settings for the current R session. For example, the
following R command specifies the use of Hartung-Knapp and
Paule-Mandel method, and the printing of prediction intervals for
any meta-analysis generated after execution of this command:
\itemize{
\item \code{settings.meta(method.random.ci = "HK", method.tau =
  "PM", prediction = TRUE)}
}
}

\subsection{Data sets}{
The following data sets are available in R package \bold{meta}.

\tabular{ll}{
\bold{Data set} \tab \bold{Description} \cr
\code{\link{Fleiss1993bin}} \tab Aspirin after myocardial infarction \cr
\code{\link{Fleiss1993cont}} \tab Mental health treatment on medical utilisation\cr
\code{\link{Olkin1995}} \tab Thrombolytic therapy after acute myocardial infarction \cr
\code{\link{Pagliaro1992}} \tab Prevention of first bleeding in cirrhosis \cr
\code{\link{amlodipine}} \tab Amlodipine for work capacity \cr
\code{\link{caffeine}} \tab Caffeine for daytime drowsiness (Cochrane Practice review) \cr
\code{\link{cisapride}} \tab Cisapride in non-ulcer dispepsia \cr
\code{\link{lungcancer}} \tab Smoking example \cr
\code{\link{smoking}} \tab Smoking example \cr
\code{\link{woodyplants}} \tab Elevated CO$_2$ and total biomass of woody plants
}

\bold{R} package \bold{metadat} has a large collection of meta-analysis data
sets.
}
}
\note{
Balduzzi et al. (2019) is the preferred citation in publications
for \bold{meta}. Type \code{citation("meta")} for a BibTeX entry of
this publication.

Type \code{help(package = "meta")} for a listing of all R functions
and datasets available in \bold{meta}. For example, results of
several meta-analyses can be combined with \code{\link{metabind}}
which is useful to generate a forest plot with results of several
subgroup analyses.

R package \bold{meta} imports R functions from \bold{metafor}
(Viechtbauer, 2010) to
\itemize{
\item estimate the between-study variance \eqn{\tau^2},
\item conduct meta-regression,
\item estimate three-level models,
\item estimate generalised linear mixed models.
}

To report problems and bugs
\itemize{
\item type \code{bug.report(package = "meta")} if you do not use
  RStudio,
\item send an email to Guido Schwarzer
  \email{guido.schwarzer@uniklinik-freiburg.de} if you use RStudio.
}

The development version of \bold{meta} is available on GitHub
\url{https://github.com/guido-s/meta/}.
}
\references{
Balduzzi S, Rücker G, Schwarzer G (2019):
How to perform a meta-analysis with R: a practical tutorial.
\emph{Evidence-Based Mental Health},
\bold{22}, 153--160

Biggerstaff BJ, Jackson D (2008):
The exact distribution of Cochran’s heterogeneity statistic in
one-way random effects meta-analysis.
\emph{Statistics in Medicine},
\bold{27}, 6093--110

DerSimonian R & Laird N (1986):
Meta-analysis in clinical trials.
\emph{Controlled Clinical Trials},
\bold{7}, 177--88

Doi SAR, Barendregt JJ, Khan S, Thalib L, Williams GM (2015):
Advances in the meta-analysis of heterogeneous clinical trials I:
The inverse variance heterogeneity model.
\emph{Contemporary Clinical Trials},
\bold{45}, 130--8

Evrenoglou T, White IR, Afach S, Mavridis D, Chaimani A. (2022):
Network meta-analysis of rare events using penalized likelihood regression.
\emph{Statistics in Medicine},
\bold{41}, 5203--19

Hartung J, Knapp G (2001a):
On tests of the overall treatment effect in meta-analysis with
normally distributed responses.
\emph{Statistics in Medicine},
\bold{20}, 1771--82

Hartung J, Knapp G (2001b):
A refined method for the meta-analysis of controlled clinical
trials with binary outcome.
\emph{Statistics in Medicine},
\bold{20}, 3875--89

Hedges LV & Olkin I (1985):
\emph{Statistical methods for meta-analysis}.
San Diego, CA: Academic Press

Higgins JPT & Thompson SG (2002):
Quantifying heterogeneity in a meta-analysis.
\emph{Statistics in Medicine},
\bold{21}, 1539--58

Higgins JPT, Thompson SG, Spiegelhalter DJ (2009):
A re-evaluation of random-effects meta-analysis.
\emph{Journal of the Royal Statistical Society: Series A},
\bold{172}, 137--59

Hunter JE & Schmidt FL (2015):
\emph{Methods of Meta-Analysis: Correcting Error and Bias in
Research Findings} (Third edition).
Thousand Oaks, CA: Sage

IntHout J, Ioannidis JPA, Borm GF (2014):
The Hartung-Knapp-Sidik-Jonkman method for random effects
meta-analysis is straightforward and considerably outperforms the
standard DerSimonian-Laird method.
\emph{BMC Medical Research Methodology},
\bold{14}, 25

IQWiG (2022):
General Methods: Version 6.1.
\url{https://www.iqwig.de/en/about-us/methods/methods-paper/}

Jackson D (2013):
Confidence intervals for the between-study variance in random
effects meta-analysis using generalised Cochran heterogeneity
statistics.
\emph{Research Synthesis Methods},
\bold{4}, 220--229

Jackson D, Law M, Rücker G, Schwarzer G (2017): 
The Hartung-Knapp modification for random-effects meta-analysis: A
useful refinement but are there any residual concerns?
\emph{Statistics in Medicine},
\bold{36}, 3923--34

Knapp G & Hartung J (2003):
Improved tests for a random effects meta-regression with a single
covariate.
\emph{Statistics in Medicine},
\bold{22}, 2693--710

Langan D, Higgins JPT, Jackson D, Bowden J, Veroniki AA,
Kontopantelis E, et al. (2019):
A comparison of heterogeneity variance estimators in simulated
random-effects meta-analyses.
\emph{Research Synthesis Methods},
\bold{10}, 83--98

Schwarzer G (2007):
meta: An R package for meta-analysis.
\emph{R News},
\bold{7}, 40--5

Schwarzer G, Carpenter JR and Rücker G (2015):
\emph{Meta-Analysis with R (Use-R!)}.
Springer International Publishing, Switzerland

Skipka G (2006):
The inclusion of the estimated inter-study variation into forest
plots for random effects meta-analysis - a suggestion for a
graphical representation [abstract].
\emph{XIV Cochrane Colloquium, Dublin}, 23-26.

Stijnen T, Hamza TH, Ozdemir P (2010):
Random effects meta-analysis of event outcome in the framework of
the generalized linear mixed model with applications in sparse
data.
\emph{Statistics in Medicine},
\bold{29}, 3046--67

Veroniki AA, Jackson D, Viechtbauer W, Bender R, Bowden J, Knapp G,
et al. (2016):
Methods to estimate the between-study variance and its uncertainty
in meta-analysis.
\emph{Research Synthesis Methods},
\bold{7}, 55--79 

Veroniki AA, Jackson D, Bender R, Kuss O, Higgins JPT, Knapp G, Salanti G
(2019):
Methods to calculate uncertainty in the estimated overall effect size from a
random-effects meta-analysis.
\emph{Research Synthesis Methods},
\bold{10}, 23--43

Van den Noortgate W, López-López JA, Marín-Martínez F, Sánchez-Meca J (2013):
Three-level meta-analysis of dependent effect sizes.
\emph{Behavior Research Methods},
\bold{45}, 576--94

Viechtbauer W (2005):
Bias and efficiency of meta-analytic variance estimators in the
random-effects model.
\emph{Journal of Educational and Behavioral Statistics},
\bold{30}, 261--93

Viechtbauer W (2007):
Confidence intervals for the amount of heterogeneity in
meta-analysis.
\emph{Statistics in Medicine},
\bold{26}, 37--52

Viechtbauer W (2010):
Conducting Meta-Analyses in R with the metafor Package.
\emph{Journal of Statistical Software},
\bold{36}, 1--48

Wiksten A, Rücker G, Schwarzer G (2016):
Hartung-Knapp method is not always conservative compared with
fixed-effect meta-analysis.
\emph{Statistics in Medicine},
\bold{35}, 2503--15
}
\seealso{
\code{\link{meta-object}}, \code{\link{meta-sm}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
\keyword{package}
