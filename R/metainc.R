#' Meta-analysis of incidence rates
#' 
#' @description
#' Calculation of fixed effect and random effects estimates (incidence
#' rate ratio or incidence rate difference) for meta-analyses with
#' event counts.  Mantel-Haenszel, Cochran, inverse variance method,
#' and generalised linear mixed model (GLMM) are available for
#' pooling. For GLMMs, the \code{\link[metafor]{rma.glmm}} function
#' from R package \bold{metafor} (Viechtbauer 2010) is called
#' internally.
#' 
#' @param event.e Number of events in experimental group.
#' @param time.e Person time at risk in experimental group.
#' @param event.c Number of events in control group.
#' @param time.c Person time at risk in control group.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information, i.e., event.e, time.e, event.c, and time.c.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"MH"},
#'   \code{"Inverse"}, \code{"Cochran"}, or \code{"GLMM"} can be
#'   abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"IRR"} or \code{"IRD"}) is to be used for pooling of
#'   studies, see Details.
#' @param incr A numerical value which is added to each cell frequency
#'   for studies with a zero cell count, see Details.
#' @param allincr A logical indicating if \code{incr} is added to each
#'   cell frequency of all studies if at least one study has a zero
#'   cell count. If FALSE (default), \code{incr} is added only to each
#'   cell frequency of studies with a zero cell count.
#' @param addincr A logical indicating if \code{incr} is added to each
#'   cell frequency of all studies irrespective of zero cell counts.
#' @param model.glmm A character string indicating which GLMM should
#'   be used. One of \code{"UM.FS"}, \code{"UM.RS"}, and
#'   \code{"CM.EL"}, see Details.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.comb The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param comb.fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param comb.random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param overall A logical indicating whether overall summaries
#'   should be reported. This argument is useful in a meta-analysis
#'   with subgroups if overall results should not be reported.
#' @param overall.hetstat A logical value indicating whether to print
#'   heterogeneity measures for overall treatment comparisons. This
#'   argument is useful in a meta-analysis with subgroups if
#'   heterogeneity statistics should only be printed on subgroup
#'   level.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param hakn A logical indicating whether the method by Hartung and
#'   Knapp should be used to adjust test statistics and confidence
#'   intervals.
#' @param adhoc.hakn A character string indicating whether an \emph{ad
#'   hoc} variance correction should be applied in the case of an
#'   arbitrarily small Hartung-Knapp variance estimate, see Details.
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau}. Either \code{"DL"}, \code{"PM"},
#'   \code{"REML"}, \code{"ML"}, \code{"HS"}, \code{"SJ"},
#'   \code{"HE"}, or \code{"EB"}, can be abbreviated.
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau}. Either \code{"QP"}, \code{"BJ"}, or \code{"J"}, or
#'   \code{""}, can be abbreviated.
#' @param tau.preset Prespecified value for the square root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance \eqn{\tau^2}.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param method.bias A character string indicating which test for
#'   funnel plot asymmetry is to be used. Either \code{"linreg"} or
#'   \code{"rank"}, can be abbreviated. See function
#'   \code{\link{metabias}}
#' @param n.e Number of observations in experimental group (optional).
#' @param n.c Number of observations in control group (optional).
#' @param backtransf A logical indicating whether results for
#'   incidence rate ratio (\code{sm = "IRR"}) should be back
#'   transformed in printouts and plots.  If TRUE (default), results
#'   will be presented as incidence rate ratios; otherwise log
#'   incidence rate ratios will be shown.
#' @param irscale A numeric defining a scaling factor for printing of
#'   incidence rate differences.
#' @param irunit A character string specifying the time unit used to
#'   calculate rates, e.g. person-years.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param label.e Label for experimental group.
#' @param label.c Label for control group.
#' @param label.left Graph label on left side of forest plot.
#' @param label.right Graph label on right side of forest plot.
#' @param byvar An optional vector containing grouping information
#'   (must be of same length as \code{event.e}).
#' @param bylab A character string with a label for the grouping
#'   variable.
#' @param print.byvar A logical indicating whether the name of the
#'   grouping variable should be printed in front of the group labels.
#' @param byseparator A character string defining the separator
#'   between label and levels of grouping variable.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if \code{incr} is added to studies with zero cell
#'   frequencies).
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}} or
#'   \code{\link[metafor]{rma.glmm}}, respectively.
#' @param \dots Additional arguments passed on to
#'   \code{\link[metafor]{rma.glmm}} function.
#' 
#' @details
#' Calculation of fixed and random effects estimates for meta-analyses
#' comparing two incidence rates.
#' 
#' The following measures of treatment effect are available:
#' \itemize{
#' \item Incidence Rate Ratio (\code{sm = "IRR"})
#' \item Incidence Rate Difference (\code{sm = "IRD"})
#' }
#' 
#' Default settings are utilised for several arguments (assignments
#' using \code{\link{gs}} function). These defaults can be changed for
#' the current R session using the \code{\link{settings.meta}}
#' function.
#' 
#' Furthermore, R function \code{\link{update.meta}} can be used to
#' rerun a meta-analysis with different settings.
#' 
#' \subsection{Meta-analysis method}{
#' 
#' By default, both fixed effect and random effects models are
#' considered (see arguments \code{comb.fixed} and
#' \code{comb.random}). If \code{method} is \code{"MH"} (default), the
#' Mantel-Haenszel method is used to calculate the fixed effect
#' estimate (Greenland & Robbins, 1985); if \code{method} is
#' \code{"Inverse"}, inverse variance weighting is used for pooling;
#' if \code{method} is \code{"Cochran"}, the Cochran method is used
#' for pooling (Bayne-Jones, 1964, Chapter 8).
#' 
#' A distinctive and frequently overlooked advantage of incidence
#' rates is that individual patient data (IPD) can be extracted from
#' count data. Accordingly, statistical methods for IPD, i.e.,
#' generalised linear mixed models, can be utilised in a meta-analysis
#' of incidence rate ratios (Stijnen et al., 2010). These methods are
#' available (argument \code{method = "GLMM"}) by calling the
#' \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor} internally.
#'
#' Three different GLMMs are available for meta-analysis of incidence
#' rate ratios using argument \code{model.glmm} (which corresponds to
#' argument \code{model} in the \code{\link[metafor]{rma.glmm}}
#' function):
#' \tabular{cl}{
#' 1. \tab Poisson regression model with fixed study effects (default)
#'  \cr
#'  \tab (\code{model.glmm = "UM.FS"}, i.e., \bold{U}nconditional
#'  \bold{M}odel - \bold{F}ixed \bold{S}tudy effects) \cr
#' 2. \tab Mixed-effects Poisson regression model with random study
#'  effects \cr
#'  \tab (\code{model.glmm = "UM.RS"}, i.e., \bold{U}nconditional
#'  \bold{M}odel - \bold{R}andom \bold{S}tudy effects) \cr
#' 3. \tab Generalised linear mixed model (conditional Poisson-Normal)
#'  \cr
#'  \tab (\code{model.glmm = "CM.EL"}, i.e., \bold{C}onditional
#'   \bold{M}odel - \bold{E}xact \bold{L}ikelihood)
#' }
#'
#' Details on these three GLMMs as well as additional arguments which
#' can be provided using argument '\code{\dots}' in \code{metainc}
#' are described in \code{\link[metafor]{rma.glmm}} where you can also
#' find information on the iterative algorithms used for estimation.
#' Note, regardless of which value is used for argument
#' \code{model.glmm}, results for two different GLMMs are calculated:
#' fixed effect model (with fixed treatment effect) and random effects
#' model (with random treatment effects).
#' }
#' 
#' \subsection{Continuity correction}{
#' 
#' For studies with a zero cell count, by default, 0.5 is added to all
#' cell frequencies of these studies (argument \code{incr}). This
#' continuity correction is used both to calculate individual study
#' results with confidence limits and to conduct meta-analysis based
#' on the inverse variance method. For Mantel-Haenszel method, Cochran
#' method, and GLMMs, nothing is added to zero cell counts.
#' Accordingly, estimates for these methods are not defined if the
#' number of events is zero in all studies either in the experimental
#' or control group.
#' }
#' 
#' \subsection{Estimation of between-study variance}{
#' 
#' The following methods to estimate the between-study variance
#' \eqn{\tau^2} are available for the inverse variance method:
#' \itemize{
#' \item DerSimonian-Laird estimator (\code{method.tau = "DL"})
#' \item Paule-Mandel estimator (\code{method.tau = "PM"})
#' \item Restricted maximum-likelihood estimator (\code{method.tau =
#'   "REML"})
#' \item Maximum-likelihood estimator (\code{method.tau = "ML"})
#' \item Hunter-Schmidt estimator (\code{method.tau = "HS"})
#' \item Sidik-Jonkman estimator (\code{method.tau = "SJ"})
#' \item Hedges estimator (\code{method.tau = "HE"})
#' \item Empirical Bayes estimator (\code{method.tau = "EB"})
#' }
#' See \code{\link{metagen}} for more information on these
#' estimators. Note, the maximum-likelihood method is utilized for
#' GLMMs.
#' }
#' 
#' \subsection{Confidence interval for the between-study variance}{
#'
#' The following methods to calculate a confidence interval for
#' \eqn{\tau^2} and \eqn{\tau} are available.
#' \tabular{ll}{
#' \bold{Argument}\tab \bold{Method} \cr 
#' \code{method.tau.ci = "J"}\tab Method by Jackson \cr
#' \code{method.tau.ci = "BJ"}\tab Method by Biggerstaff and Jackson \cr
#' \code{method.tau.ci = "QP"}\tab Q-Profile method
#' }
#' See \code{\link{metagen}} for more information on these
#' methods. For GLMMs, no confidence intervals for \eqn{\tau^2} and
#' \eqn{\tau} are calculated. Likewise, no confidence intervals for
#' \eqn{\tau^2} and \eqn{\tau} are calculated if \code{method.tau.ci =
#' ""}. 
#' }
#' 
#' \subsection{Hartung-Knapp method}{
#' 
#' Hartung and Knapp (2001a,b) proposed an alternative method for
#' random effects meta-analysis based on a refined variance estimator
#' for the treatment estimate. Simulation studies (Hartung and Knapp,
#' 2001a,b; IntHout et al., 2014; Langan et al., 2019) show improved
#' coverage probabilities compared to the classic random effects
#' method.
#'
#' In rare settings with very homogeneous treatment estimates, the
#' Hartung-Knapp variance estimate can be arbitrarily small resulting
#' in a very narrow confidence interval (Knapp and Hartung, 2003;
#' Wiksten et al., 2016). In such cases, an \emph{ad hoc} variance
#' correction has been proposed by utilising the variance estimate
#' from the classic random effects model (Knapp and Hartung,
#' 2003). Argument \code{adhoc.hakn} can be used to choose the
#' \emph{ad hoc} method:
#' \tabular{ll}{
#' \bold{Argument}\tab \bold{\emph{Ad hoc} method} \cr 
#' \code{adhoc.hakn = ""}\tab not used \cr
#' \code{adhoc.hakn = "se"}\tab used if HK standard error is smaller than
#'  standard error \cr
#'  \tab from classic random effects model (Knapp and Hartung, 2003) \cr
#' \code{adhoc.hakn = "ci"}\tab used if HK confidence interval is
#'  narrower than CI from \cr
#'  \tab classic random effects model with DL estimator (IQWiG, 2020)
#' }
#' }
#' 
#' \subsection{Prediction interval}{
#' 
#' A prediction interval for the proportion in a new study (Higgins et
#' al., 2009) is calculated if arguments \code{prediction} and
#' \code{comb.random} are \code{TRUE}. Note, the definition of
#' prediction intervals varies in the literature. This function
#' implements equation (12) of Higgins et al., (2009) which proposed a
#' \emph{t} distribution with \emph{K-2} degrees of freedom where
#' \emph{K} corresponds to the number of studies in the meta-analysis.
#' 
#' For GLMMs, a method similar to Knapp and Hartung (2003) is
#' implemented, see description of argument \code{tdist} in
#' \code{\link[metafor]{rma.glmm}}.
#' }
#'
#' \subsection{Subgroup analysis}{
#' 
#' Argument \code{byvar} can be used to conduct subgroup analysis for
#' a categorical covariate. The \code{\link{metareg}} function can be
#' used instead for more than one categorical covariate or continuous
#' covariates.
#' }
#' 
#' \subsection{Exclusion of studies from meta-analysis}{
#'
#' Arguments \code{subset} and \code{exclude} can be used to exclude
#' studies from the meta-analysis. Studies are removed completely from
#' the meta-analysis using argument \code{subset}, while excluded
#' studies are shown in printouts and forest plots using argument
#' \code{exclude} (see Examples in \code{\link{metagen}}).
#' Meta-analysis results are the same for both arguments.
#' }
#' 
#' \subsection{Presentation of meta-analysis results}{
#' 
#' Internally, both fixed effect and random effects models are
#' calculated regardless of values choosen for arguments
#' \code{comb.fixed} and \code{comb.random}. Accordingly, the estimate
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"meta"} even if
#' argument \code{comb.random = FALSE}. However, all functions in R
#' package \bold{meta} will adequately consider the values for
#' \code{comb.fixed} and \code{comb.random}. E.g. function
#' \code{\link{print.meta}} will not print results for the random
#' effects model if \code{comb.random = FALSE}.
#' }
#' 
#' @return
#' An object of class \code{c("metainc", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#' \item{event.e, time.e, event.c, time.c, studlab, exclude,}{As
#'   defined above.}
#' \item{sm, method, incr, allincr, addincr, model.glmm, warn,}{As
#'   defined above.}
#' \item{level, level.comb, comb.fixed, comb.random,}{As defined
#'   above.}
#' \item{overall, overall.hetstat,}{As defined above.}
#' \item{hakn, adhoc.hakn, method.tau, method.tau.ci,}{As defined above.}
#' \item{tau.preset, TE.tau, method.bias,}{As defined above.}
#' \item{tau.common, title, complab, outclab,}{As defined above.}
#' \item{label.e, label.c, label.left, label.right,}{As defined
#'   above.}
#' \item{byvar, bylab, print.byvar, byseparator}{As defined above.}
#' \item{TE, seTE}{Estimated treatment effect and standard error of
#'   individual studies.}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies.}
#' \item{zval, pval}{z-value and p-value for test of treatment effect
#'   for individual studies.}
#' \item{w.fixed, w.random}{Weight of individual studies (in fixed and
#'   random effects model).}
#' \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
#'   standard error (fixed effect model).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (fixed effect model).}
#' \item{zval.fixed, pval.fixed}{z-value and p-value for test of
#'   overall treatment effect (fixed effect model).}
#' \item{TE.random, seTE.random}{Estimated overall treatment effect
#'   and standard error (random effects model).}
#' \item{lower.random, upper.random}{Lower and upper confidence
#'   interval limits (random effects model).}
#' \item{zval.random, pval.random}{z-value or t-value and
#'   corresponding p-value for test of overall treatment effect
#'   (random effects model).}
#' \item{prediction, level.predict}{As defined above.}
#' \item{seTE.predict}{Standard error utilised for prediction
#'   interval.}
#' \item{lower.predict, upper.predict}{Lower and upper limits of
#'   prediction interval.}
#' \item{k}{Number of studies combined in meta-analysis.}
#' \item{Q}{Heterogeneity statistic Q.}
#' \item{df.Q}{Degrees of freedom for heterogeneity statistic.}
#' \item{pval.Q}{P-value of heterogeneity test.}
#' \item{Q.LRT}{Heterogeneity statistic for likelihood-ratio test
#'   (only if \code{method = "GLMM"}).}
#' \item{df.Q.LRT}{Degrees of freedom for likelihood-ratio test}
#' \item{pval.Q.LRT}{P-value of likelihood-ratio test.}
#' \item{tau2}{Between-study variance \eqn{\tau^2}.}
#' \item{se.tau2}{Standard error of \eqn{\tau^2}.}
#' \item{lower.tau2, upper.tau2}{Lower and upper limit of confidence
#'   interval for \eqn{\tau^2}.}
#' \item{tau}{Square-root of between-study variance \eqn{\tau}.}
#' \item{lower.tau, upper.tau}{Lower and upper limit of confidence
#'   interval for \eqn{\tau}.}
#' \item{H}{Heterogeneity statistic H.}
#' \item{lower.H, upper.H}{Lower and upper confidence limit for
#'  heterogeneity statistic H.}
#' \item{I2}{Heterogeneity statistic I\eqn{^2}.}
#' \item{lower.I2, upper.I2}{Lower and upper confidence limit for
#'   heterogeneity statistic I\eqn{^2}.}
#' \item{Rb}{Heterogeneity statistic R\eqn{_b}.}
#' \item{lower.Rb, upper.Rb}{Lower and upper confidence limit for
#'   heterogeneity statistic R\eqn{_b}.}
#' \item{sparse}{Logical flag indicating if any study included in
#'   meta-analysis has any zero cell frequencies.}
#' \item{incr.event}{Increment added to number of events.}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn = TRUE}).}
#' \item{k.MH}{Number of studies combined in meta-analysis using
#'   Mantel-Haenszel method.}
#' \item{bylevs}{Levels of grouping variable - if \code{byvar} is not
#'   missing.}
#' \item{TE.fixed.w, seTE.fixed.w}{Estimated treatment effect and
#'   standard error in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{lower.fixed.w, upper.fixed.w}{Lower and upper confidence
#'   interval limits in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{zval.fixed.w, pval.fixed.w}{z-value and p-value for test of
#'   treatment effect in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}
#' \item{TE.random.w, seTE.random.w}{Estimated treatment effect and
#'   standard error in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{zval.random.w, pval.random.w}{z-value or t-value and
#'   corresponding p-value for test of treatment effect in subgroups
#'   (random effects model) - if \code{byvar} is not missing.}
#' \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
#'   random effects model) - if \code{byvar} is not missing.}
#' \item{df.hakn.w}{Degrees of freedom for test of treatment effect
#'   for Hartung-Knapp method in subgroups - if \code{byvar} is not
#'   missing and \code{hakn = TRUE}.}
#' \item{event.e.w}{Number of events in experimental group in
#'   subgroups - if \code{byvar} is not missing.}
#' \item{time.e.w}{Total person time in subgroups (experimental
#'   group) - if \code{byvar} is not missing.}
#' \item{n.e.w}{Number of observations in experimental group in
#'   subgroups - if \code{byvar} is not missing.}
#' \item{event.c.w}{Number of events in control group in subgroups -
#'   if \code{byvar} is not missing.}
#' \item{time.c.w}{Total person time in subgroups (control group) - if
#'   \code{byvar} is not missing.}
#' \item{n.c.w}{Number of observations in control group in subgroups -
#'   if \code{byvar} is not missing.}
#' \item{k.w}{Number of studies combined within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{k.all.w}{Number of all studies in subgroups - if \code{byvar}
#'   is not missing.}
#' \item{Q.w.fixed}{Overall within subgroups heterogeneity statistic Q
#'   (based on fixed effect model) - if \code{byvar} is not missing.}
#' \item{Q.w.random}{Overall within subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{byvar} is not
#'   missing (only calculated if argument \code{tau.common} is TRUE).}
#' \item{df.Q.w}{Degrees of freedom for test of overall within
#'   subgroups heterogeneity - if \code{byvar} is not missing.}
#' \item{pval.Q.w.fixed}{P-value of within subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{byvar} is
#'   not missing.}
#' \item{pval.Q.w.random}{P-value of within subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{byvar} is
#'   not missing.}
#' \item{Q.b.fixed}{Overall between subgroups heterogeneity statistic
#'   Q (based on fixed effect model) - if \code{byvar} is not
#'   missing.}
#' \item{Q.b.random}{Overall between subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{byvar} is not
#'   missing.}
#' \item{df.Q.b}{Degrees of freedom for test of overall between
#'   subgroups heterogeneity - if \code{byvar} is not missing.}
#' \item{pval.Q.b.fixed}{P-value of between subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{byvar} is
#'   not missing.}
#' \item{pval.Q.b.random}{P-value of between subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{byvar} is
#'   not missing.}
#' \item{tau.w}{Square-root of between-study variance within subgroups
#'   - if \code{byvar} is not missing.}
#' \item{H.w}{Heterogeneity statistic H within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.H.w, upper.H.w}{Lower and upper confidence limit for
#'   heterogeneity statistic H within subgroups - if \code{byvar} is
#'   not missing.}
#' \item{I2.w}{Heterogeneity statistic I\eqn{^2} within subgroups - if
#'   \code{byvar} is not missing.}
#' \item{lower.I2.w, upper.I2.w}{Lower and upper confidence limit for
#'   heterogeneity statistic I\eqn{^2} within subgroups - if \code{byvar} is
#'   not missing.}
#' \item{keepdata}{As defined above.}
#' \item{data}{Original data (set) used in function call (if
#'   \code{keepdata = TRUE}).}
#' \item{subset}{Information on subset of original data used in
#'   meta-analysis (if \code{keepdata = TRUE}).}
#' \item{.glmm.fixed}{GLMM object generated by call of
#'   \code{\link[metafor]{rma.glmm}} function (fixed effect model).}
#' \item{.glmm.random}{GLMM object generated by call of
#'   \code{\link[metafor]{rma.glmm}} function (random effects model).}
#' \item{call}{Function call.}
#' \item{version}{Version of R package \bold{meta} used to create
#'   object.}
#' \item{version.metafor}{Version of R package \bold{metafor} used for
#'   GLMMs.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{metabin}}, \code{\link{update.meta}},
#'   \code{\link{print.meta}}
#' 
#' @references
#' Bayne-Jones S et al. (1964):
#' Smoking and Health: Report of the Advisory Committee to the Surgeon
#' General of the United States.
#' U-23 Department of Health, Education, and Welfare.
#' Public Health Service Publication No. 1103.
#' 
#' DerSimonian R & Laird N (1986):
#' Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials},
#' \bold{7}, 177--88
#' 
#' Greenland S & Robins JM (1985):
#' Estimation of a common effect parameter from sparse follow-up data.
#' \emph{Biometrics},
#' \bold{41}, 55--68
#' 
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009):
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{172}, 137--59
#' 
#' Knapp G & Hartung J (2003):
#' Improved tests for a random effects meta-regression with a single
#' covariate.
#' \emph{Statistics in Medicine},
#' \bold{22}, 2693--710
#' 
#' Paule RC & Mandel J (1982):
#' Consensus values and weighting factors.
#'\emph{Journal of Research of the National Bureau of Standards},
#' \bold{87}, 377--85
#' 
#' Stijnen T, Hamza TH, Ozdemir P (2010):
#' Random effects meta-analysis of event outcome in the framework of
#' the generalized linear mixed model with applications in sparse
#' data.
#' \emph{Statistics in Medicine},
#' \bold{29}, 3046--67
#' 
#' Viechtbauer W (2010):
#' Conducting Meta-Analyses in R with the Metafor Package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#' 
#' Wiksten A, Rücker G, Schwarzer G (2016):
#' Hartung-Knapp method is not always conservative compared with
#' fixed-effect meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{35}, 2503--15
#' 
#' @examples
#' data(smoking)
#' m1 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
#'               data = smoking, studlab = study)
#' print(m1, digits = 2)
#' 
#' m2 <- update(m1, method = "Cochran")
#' print(m2, digits = 2)
#' 
#' data(lungcancer)
#' m3 <- metainc(d.smokers, py.smokers,
#'               d.nonsmokers, py.nonsmokers,
#'               data = lungcancer, studlab = study)
#' print(m3, digits = 2)
#' 
#' # Redo Cochran meta-analysis with inflated standard errors
#' #
#' # All cause mortality
#' #
#' TEa <- log((smoking$d.smokers/smoking$py.smokers) /
#'            (smoking$d.nonsmokers/smoking$py.nonsmokers))
#' seTEa <- sqrt(1 / smoking$d.smokers + 1 / smoking$d.nonsmokers +
#'               2.5 / smoking$d.nonsmokers)
#' metagen(TEa, seTEa, sm = "IRR", studlab = smoking$study)
#' 
#' # Lung cancer mortality
#' #
#' TEl <- log((lungcancer$d.smokers/lungcancer$py.smokers) /
#'            (lungcancer$d.nonsmokers/lungcancer$py.nonsmokers))
#' seTEl <- sqrt(1 / lungcancer$d.smokers + 1 / lungcancer$d.nonsmokers +
#'               2.25 / lungcancer$d.nonsmokers)
#' metagen(TEl, seTEl, sm = "IRR", studlab = lungcancer$study)
#'
#' \dontrun{
#' # Meta-analysis using generalised linear mixed models
#' # (only if R packages 'metafor' and 'lme4' are available)
#' 
#' # Poisson regression model (fixed study effects)
#' #
#' m4 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
#'               data = smoking, studlab = study, method = "GLMM")
#' m4
#' 
#' # Mixed-effects Poisson regression model (random study effects)
#' #
#' update(m4, model.glmm = "UM.RS", nAGQ = 1)
#' #
#' # Generalised linear mixed model (conditional Poisson-Normal)
#' #
#' update(m4, model.glmm = "CM.EL")
#' }
#' 
#' @export metainc


metainc <- function(event.e, time.e, event.c, time.c, studlab,
                    ##
                    data = NULL, subset = NULL, exclude = NULL,
                    method = "MH",
                    ##
                    sm = gs("sminc"),
                    ##
                    incr = gs("incr"), allincr = gs("allincr"),
                    addincr = gs("addincr"),
                    model.glmm = "UM.FS",
                    ##
                    level = gs("level"), level.comb = gs("level.comb"),
                    comb.fixed = gs("comb.fixed"),
                    comb.random = gs("comb.random"),
                    overall = comb.fixed | comb.random,
                    overall.hetstat = comb.fixed | comb.random,
                    ##
                    hakn = gs("hakn"), adhoc.hakn = gs("adhoc.hakn"),
                    method.tau =
                      ifelse(!is.na(charmatch(tolower(method), "glmm",
                                              nomatch = NA)),
                             "ML", gs("method.tau")),
                    method.tau.ci = if (method.tau == "DL") "J" else "QP",
                    tau.preset = NULL, TE.tau = NULL,
                    tau.common = gs("tau.common"),
                    ##
                    prediction = gs("prediction"),
                    level.predict = gs("level.predict"),
                    ##
                    method.bias = gs("method.bias"),
                    ##
                    n.e = NULL, n.c = NULL,
                    ##
                    backtransf = gs("backtransf"),
                    irscale = 1, irunit = "person-years",
                    title = gs("title"), complab = gs("complab"),
                    outclab = "",
                    label.e = gs("label.e"), label.c = gs("label.c"),
                    label.left = gs("label.left"),
                    label.right = gs("label.right"),
                    ##
                    byvar, bylab, print.byvar = gs("print.byvar"),
                    byseparator = gs("byseparator"),
                    ##
                    keepdata = gs("keepdata"),
                    warn = gs("warn"),
                    ##
                    control = NULL,
                    ...
                    ) {
  
  
  ##
  ##
  ## (1) Check arguments
  ##
  ##
  chknull(sm)
  sm <- setchar(sm, c("IRR", "IRD"))
  ##
  chklevel(level)
  chklevel(level.comb)
  chklogical(comb.fixed)
  chklogical(comb.random)
  chklogical(overall)
  chklogical(overall.hetstat)
  ##
  chklogical(hakn)
  adhoc.hakn <- setchar(adhoc.hakn, .settings$adhoc4hakn)
  method.tau <- setchar(method.tau, .settings$meth4tau)
  method.tau.ci <- setchar(method.tau.ci, .settings$meth4tau.ci)
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  method.bias <- setchar(method.bias, .settings$meth4bias)
  ##
  chklogical(backtransf)
  ##
  chknumeric(irscale, length = 1)
  chkchar(irunit)
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metainc"
  ##
  if (sm != "IRD" & irscale != 1) {
    warning("Argument 'irscale' only considered for incidence rate differences.")
    irscale <- 1
  }
  ##
  method <- setchar(method, c("Inverse", "MH", "Cochran", "GLMM"))
  is.glmm <- method == "GLMM"
  ##
  chklogical(allincr)
  chklogical(addincr)
  model.glmm <- setchar(model.glmm, c("UM.FS", "UM.RS", "CM.EL"))
  chklogical(warn)
  ##
  if (is.glmm & sm != "IRR")
    stop("Generalised linear mixed models only possible with argument 'sm = \"IRR\"'.")
  ##
  if (is.glmm & method.tau != "ML")
    stop("Generalised linear mixed models only possible with argument 'method.tau = \"ML\"'.")
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch 'event.e', 'time.e', 'event.c', 'time.c', 'n.e', and 'n.c' from data:
  ##
  event.e <- eval(mf[[match("event.e", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.e)
  k.All <- length(event.e)
  ##
  time.e <- eval(mf[[match("time.e", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(time.e)
  ##
  event.c <- eval(mf[[match("event.c", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.c)
  ##
  time.c <- eval(mf[[match("time.c", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknull(time.c)
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  ##
  ## Catch 'incr' from data:
  ##
  if (!missing(incr))
    incr <- eval(mf[[match("incr", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknumeric(incr, min = 0)
  ##
  ## Catch 'studlab', 'byvar', 'subset', and 'exclude' from data:
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  ##
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  by <- !is.null(byvar)
  ##
  subset <- eval(mf[[match("subset", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  missing.subset <- is.null(subset)
  ##
  exclude <- eval(mf[[match("exclude", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  missing.exclude <- is.null(exclude)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  chklength(time.e, k.All, fun)
  chklength(event.c, k.All, fun)
  chklength(time.c, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (length(incr) > 1)
    chklength(incr, k.All, fun)
  ##
  if (by)
    chklength(byvar, k.All, fun)
  ##
  if (!is.null(n.e))
    chklength(n.e, k.All, fun)
  if (!is.null(n.c))
    chklength(n.c, k.All, fun)
  ##
  ## Additional checks
  ##
  if (is.glmm) {
    if (!is.null(TE.tau)) {
      if (warn)
        warning("Argument 'TE.tau' not considered for GLMM.")
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)) {
      if (warn)
        warning("Argument 'tau.preset' not considered for GLMM.")
      tau.preset <- NULL
    }
  }
  if (!by & tau.common) {
    warning("Value for argument 'tau.common' set to FALSE as argument 'byvar' is missing.")
    tau.common <- FALSE
  }
  if (by & !tau.common & !is.null(tau.preset)) {
    warning("Argument 'tau.common' set to TRUE as argument tau.preset is not NULL.")
    tau.common <- TRUE
  }
  
  
  ##
  ##
  ## (4) Subset, exclude studies, and subgroups
  ##
  ##
  if (!missing.subset)
    if ((is.logical(subset) & (sum(subset) > k.All)) ||
        (length(subset) > k.All))
      stop("Length of subset is larger than number of studies.")
  ##
  if (!missing.exclude) {
    if ((is.logical(exclude) & (sum(exclude) > k.All)) ||
        (length(exclude) > k.All))
      stop("Length of argument 'exclude' is larger than number of studies.")
    ##
    exclude2 <- rep(FALSE, k.All)
    exclude2[exclude] <- TRUE
    exclude <- exclude2
  }
  else
    exclude <- rep(FALSE, k.All)
  ##
  if (by) {
    chkmiss(byvar)
    byvar.name <- byvarname(mf[[match("byvar", names(mf))]])
    bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  
  
  ##
  ##
  ## (5) Store complete dataset in list object data
  ##     (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata) {
    if (nulldata)
      data <- data.frame(.event.e = event.e)
    else
      data$.event.e <- event.e
    ##
    data$.time.e <- time.e
    data$.event.c <- event.c
    data$.time.c <- time.c
    data$.studlab <- studlab
    ##
    data$.incr <- incr
    ##
    if (by)
      data$.byvar <- byvar
    ##
    if (!missing.subset) {
      if (length(subset) == dim(data)[1])
        data$.subset <- subset
      else {
        data$.subset <- FALSE
        data$.subset[subset] <- TRUE
      }
    }
    ##
    if (!missing.exclude)
      data$.exclude <- exclude
    ##
    if (!is.null(n.e))
      data$.n.e <- n.e
    if (!is.null(n.e))
      data$.n.c <- n.c
  }  
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    event.e <- event.e[subset]
    time.e <- time.e[subset]
    event.c <- event.c[subset]
    time.c <- time.c[subset]
    studlab <- studlab[subset]
    ##
    exclude <- exclude[subset]
    ##
    if (length(incr) > 1)
      incr <- incr[subset]
    ##
    if (by)
      byvar <- byvar[subset]
    ##
    if (!is.null(n.e))
      n.e <- n.e[subset]
    if (!is.null(n.c))
      n.c <- n.c[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(event.e)
  ##
  if (k.all == 0)
    stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1) {
    comb.fixed  <- FALSE
    comb.random <- FALSE
    prediction  <- FALSE
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(event.e, 0)
  chknumeric(time.e, 0, zero = TRUE)
  chknumeric(event.c, 0)
  chknumeric(time.c, zero = TRUE)
  ##
  ## Recode integer as numeric:
  ##
  event.e <- int2num(event.e)
  time.e  <- int2num(time.e)
  event.c <- int2num(event.c)
  time.c  <- int2num(time.c)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  sel <- switch(sm,
                IRD = event.e == 0 | event.c == 0,
                IRR = event.e == 0 | event.c == 0)
  ##
  ## Sparse computation
  ##
  sparse <- any(sel, na.rm = TRUE)
  ##
  if (is.glmm & sparse)
    if ((!missing(incr) & any(incr != 0)) |
        (!missing(allincr) & allincr ) |
        (!missing(addincr) & addincr)
        )
      warning("Note, for method = \"GLMM\", continuity correction only used to calculate individual study results.")
  ##
  if (addincr)
    incr.event <- if (length(incr) == 1) rep(incr, k.all) else incr
  else
    if (sparse)
      if (allincr)
        incr.event <- if (length(incr) == 1) rep(incr, k.all) else incr
      else
        incr.event <- incr * sel
    else
      incr.event <- rep(0, k.all)
  ##  
  if (sm == "IRR") {
    TE <- log(((event.e + incr.event) / time.e) / ((event.c + incr.event) / time.c))
    seTE <- sqrt(1 / (event.e + incr.event) + 1 / (event.c + incr.event))
  }
  else if (sm == "IRD") {
    TE <- event.e / time.e - event.c / time.c
    seTE <- sqrt((event.e + incr.event) / time.e^2 + (event.c + incr.event) / time.c^2)
  }
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(event.e[!exclude]) & !is.na(event.c[!exclude]) &
           !is.na(time.e[!exclude]) & !is.na(time.c[!exclude]))
  ##
  if (method == "MH") {
    ##
    ## Greenland, Robins (1985)
    ## 
    x.k <- event.e
    y.k <- event.c
    n.k <- time.e
    m.k <- time.c
    ##
    N.k <- n.k + m.k
    t.k <- x.k + y.k
    ##
    if (sm == "IRR") {
      D <- n.k * m.k * t.k / N.k^2
      R <- x.k * m.k / N.k
      S <- y.k * n.k / N.k
      ##
      D[exclude] <- R[exclude] <- S[exclude] <- 0
      ##
      w.fixed <- S
      TE.fixed <- log(sum(R, na.rm = TRUE) / sum(S, na.rm = TRUE))
      seTE.fixed <- sqrt(sum(D, na.rm = TRUE) / (sum(R, na.rm = TRUE) *
                                                   sum(S, na.rm = TRUE)))
    }
    else if (sm == "IRD") {
      L <- (x.k * m.k^2 + y.k * n.k^2) / N.k^2
      S <- n.k * m.k / N.k
      ##
      L[exclude] <- S[exclude] <- 0
      ##
      w.fixed <- S
      TE.fixed <- weighted.mean(TE, w.fixed, na.rm = TRUE)
      seTE.fixed <- sqrt(sum(L, na.rm = TRUE) / sum(S, na.rm = TRUE)^2)
    }
  }
  ##
  else if (method == "Cochran") {
    ##
    ## Smoking and Health - Report of the Advisory Committee to the
    ## Surgeon General of the Public Health Service,
    ## Chapter 8
    ## 
    if (sm == "IRR") {
      w.fixed <- event.c * time.e / time.c
      w.fixed[exclude] <- 0
      TE.fixed <- weighted.mean(TE, w.fixed)
      seTE.fixed <- sqrt(1 / sum(event.e) + 1 / sum(event.c))
    }
    else if (sm == "IRD") {
      warning("Cochran method only available for Incidence Rate Ratio (sm = \"IRR\")")
      return(NULL)
    }
  }
  else if (is.glmm) {
    glmm.fixed <- rma.glmm(x1i = event.e[!exclude], t1i = time.e[!exclude],
                           x2i = event.c[!exclude], t2i = time.c[!exclude],
                           method = "FE", test = ifelse(hakn, "t", "z"),
                           level = 100 * level.comb,
                           measure = "IRR", model = model.glmm,
                           control = control,
                           ...)
    ##
    TE.fixed   <- as.numeric(glmm.fixed$b)
    seTE.fixed <- as.numeric(glmm.fixed$se)
    ##
    w.fixed <- rep(NA, length(time.e))
  }
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
               ##
               sm = sm,
               level = level,
               level.comb = level.comb,
               comb.fixed = comb.fixed,
               comb.random = comb.random,
               overall = overall,
               overall.hetstat = overall.hetstat,
               ##
               hakn = hakn, adhoc.hakn = adhoc.hakn,
               method.tau = method.tau, method.tau.ci = method.tau.ci,
               tau.preset = tau.preset,
               TE.tau = if (method == "Inverse") TE.tau else TE.fixed,
               tau.common = FALSE,
               ##
               prediction = prediction,
               level.predict = level.predict,
               ##
               method.bias = method.bias,
               ##
               backtransf = backtransf,
               title = title, complab = complab, outclab = outclab,
               label.e = label.e, label.c = label.c,
               label.left = label.left, label.right = label.right,
               ##
               keepdata = FALSE,
               warn = warn,
               ##
               control = control)
  ##
  if (by & tau.common) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau, "",
                   if (method == "Inverse") TE.tau else TE.fixed,
                   level.comb, byvar, control)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(event.e = event.e, time.e = time.e,
              event.c = event.c, time.c = time.c,
              method = method,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              sparse = sparse,
              allincr = allincr, addincr = addincr,
              incr.event = incr.event,
              k.MH = if (method == "MH") sum(w.fixed > 0) else NA)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  m$method <- NULL
  ##
  res <- c(res, m)
  ##
  ## Add data
  ##
  res$n.e <- n.e
  res$n.c <- n.c
  res$TE.tau <- TE.tau
  ##
  res$irscale <- irscale
  res$irunit  <- irunit
  ##
  res$call <- match.call()
  ##
  if (method %in% c("MH", "Cochran", "GLMM")) {
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level = level.comb)
    ##
    res$TE.fixed <- TE.fixed
    res$seTE.fixed <- seTE.fixed
    res$w.fixed <- w.fixed
    res$lower.fixed <- ci.f$lower
    res$upper.fixed <- ci.f$upper
    res$statistic.fixed <- ci.f$statistic
    res$pval.fixed <- ci.f$p
  }
  ##
  if (is.glmm) {
    ##
    if (sum(!exclude) > 1)
      glmm.random <- rma.glmm(x1i = event.e[!exclude], t1i = time.e[!exclude],
                              x2i = event.c[!exclude], t2i = time.c[!exclude],
                              method = method.tau,
                              test = ifelse(hakn, "t", "z"),
                              level = 100 * level.comb,
                              measure = "IRR", model = model.glmm,
                              control = control,
                              ...)
    else {
      ##
      ## Fallback to fixed effect model due to small number of studies
      ##
      glmm.random <- glmm.fixed
    }
    ##
    TE.random   <- as.numeric(glmm.random$b)
    seTE.random <- as.numeric(glmm.random$se)
    ##
    ci.r <- ci(TE.random, seTE.random, level = level.comb)
    ##
    res$w.random <- rep(NA, length(event.e))
    ##
    res$TE.random <- TE.random
    res$seTE.random <- seTE.random
    res$lower.random <- ci.r$lower
    res$upper.random <- ci.r$upper
    res$statistic.random <- ci.r$statistic
    res$pval.random <- ci.r$p
    ##
    ## Prediction interval
    ##
    if (k >= 3) {
      tau2.calc <- if (is.na(glmm.random$tau2)) 0 else glmm.random$tau2
      seTE.predict <- sqrt(seTE.random^2 + tau2.calc)
      ci.p <- ci(TE.random, seTE.predict, level.predict, k - 2)
      res$seTE.predict <- seTE.predict
      res$lower.predict <- ci.p$lower
      res$upper.predict <- ci.p$upper
    }
    else {
      res$seTE.predict <- NA
      res$lower.predict <- NA
      res$upper.predict <- NA
    }
    ##
    res$model.glmm <- model.glmm
    ##
    res$Q <- glmm.random$QE.Wld
    res$df.Q <- glmm.random$QE.df
    res$pval.Q <- pvalQ(res$Q, res$df.Q)
    ##
    res$Q.LRT      <- glmm.random$QE.LRT
    res$df.Q.LRT   <- res$df.Q
    res$pval.Q.LRT <- pvalQ(res$Q.LRT, res$df.Q.LRT)
    ##
    if (k > 1) {
      res$tau <- sqrt(glmm.random$tau2)
      res$tau2 <- glmm.random$tau2
      res$se.tau2 <- glmm.random$se.tau2
    }
    else
      res$se.tau2 <- NA
    ##
    res$lower.tau2 <- NA
    res$upper.tau2 <- NA
    ##
    res$lower.tau <- NA
    res$upper.tau <- NA
    ##
    res$method.tau.ci <- ""
    res$sign.lower.tau <- ""
    res$sign.upper.tau <- ""
    ##
    res$H <- sqrt(glmm.random$H2)
    res$lower.H <- NA
    res$upper.H <- NA
    ##
    res$I2 <- glmm.random$I2 / 100
    res$lower.I2 <- NA
    res$upper.I2 <- NA
    ##
    res$Rb <- NA
    res$lower.Rb <- NA
    res$upper.Rb <- NA
    ##
    res$.glmm.fixed  <- glmm.fixed
    res$.glmm.random <- glmm.random
    res$version.metafor <- packageDescription("metafor")$Version
  }
  ##
  if (keepdata) {
    res$data <- data
    if (!missing.subset)
      res$subset <- subset
  }
  ##
  class(res) <- c(fun, "meta")
  ##
  ## Add results from subgroup analysis
  ##
  if (by) {
    res$byvar <- byvar
    res$bylab <- bylab
    res$print.byvar <- print.byvar
    res$byseparator <- byseparator
    res$tau.common <- tau.common
    ##
    if (!tau.common) {
      res <- c(res, subgroup(res))
      res$tau2.resid <- res$lower.tau2.resid <- res$upper.tau2.resid <- NA
      res$tau.resid <- res$lower.tau.resid <- res$upper.tau.resid <- NA
    }
    else if (!is.null(tau.preset)) {
      res <- c(res, subgroup(res, tau.preset))
      res$tau2.resid <- res$lower.tau2.resid <- res$upper.tau2.resid <- NA
      res$tau.resid <- res$lower.tau.resid <- res$upper.tau.resid <- NA
    }
    else {
      if (is.glmm) {
        res <- c(res, subgroup(res, NULL,
                               factor(res$byvar, bylevs(res$byvar)), ...))
        res$tau2.resid <- res$lower.tau2.resid <- res$upper.tau2.resid <- NA
        res$tau.resid <- res$lower.tau.resid <- res$upper.tau.resid <- NA
      }
      else {
        res <- c(res, subgroup(res, hcc$tau))
        res$Q.w.random <- hcc$Q
        res$df.Q.w.random <- hcc$df.Q
        res$tau2.resid <- hcc$tau2
        res$lower.tau2.resid <- hcc$lower.tau2
        res$upper.tau2.resid <- hcc$upper.tau2
        res$tau.resid <- hcc$tau
        res$lower.tau.resid <- hcc$lower.tau
        res$upper.tau.resid <- hcc$upper.tau
        res$sign.lower.tau.resid <- hcc$sign.lower.tau
        res$sign.upper.tau.resid <- hcc$sign.upper.tau
      }
    }
    ##
    if (!tau.common || method.tau == "DL") {
      ci.H.resid <- calcH(res$Q.w.fixed, res$df.Q.w, level.comb)
      ##
      res$H.resid <- ci.H.resid$TE
      res$lower.H.resid <- ci.H.resid$lower
      res$upper.H.resid <- ci.H.resid$upper
      ##
      ci.I2.resid <- isquared(res$Q.w.fixed, res$df.Q.w, level.comb)
      ##
      res$I2.resid <- ci.I2.resid$TE
      res$lower.I2.resid <- ci.I2.resid$lower
      res$upper.I2.resid <- ci.I2.resid$upper
    }
    else {
      if (is.glmm) {
        ci.H.resid <- calcH(res$Q.w.fixed, res$df.Q.w, level.comb)
        ##
        res$H.resid <- ci.H.resid$TE
        res$lower.H.resid <- ci.H.resid$lower
        res$upper.H.resid <- ci.H.resid$upper
        ##
        ci.I2.resid <- isquared(res$Q.w.fixed, res$df.Q.w, level.comb)
        ##
        res$I2.resid <- ci.I2.resid$TE
        res$lower.I2.resid <- ci.I2.resid$lower
        res$upper.I2.resid <- ci.I2.resid$upper
      }
      else {
        res$H.resid <- hcc$H.resid
        res$lower.H.resid <- hcc$lower.H.resid
        res$upper.H.resid <- hcc$upper.H.resid
        ##
        res$I2.resid <- hcc$I2.resid
        res$lower.I2.resid <- hcc$lower.I2.resid
        res$upper.I2.resid <- hcc$upper.I2.resid
      }
    }
    ##
    res$event.w <- NULL
    res$n.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
