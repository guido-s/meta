#' Meta-analysis of single proportions
#' 
#' @description
#' Calculation of an overall proportion from studies reporting a
#' single proportion. Inverse variance method and generalised linear
#' mixed model (GLMM) are available for pooling. For GLMMs, the
#' \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor} (Viechtbauer 2010) is called internally.
#' 
#' @param event Number of events.
#' @param n Number of observations.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information, i.e., event and n.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"Inverse"} and
#'   \code{"GLMM"}, can be abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"PLOGIT"}, \code{"PAS"}, \code{"PFT"}, \code{"PLN"}, or
#'   \code{"PRAW"}) is to be used for pooling of studies, see Details.
#' @param incr A numeric which is added to event number and sample
#'   size of studies with zero or all events, i.e., studies with an
#'   event probability of either 0 or 1.
#' @param allincr A logical indicating if \code{incr} is considered
#'   for all studies if at least one study has either zero or all
#'   events. If FALSE (default), \code{incr} is considered only in
#'   studies with zero or all events.
#' @param addincr A logical indicating if \code{incr} is used for all
#'   studies irrespective of number of events.
#' @param method.ci A character string indicating which method is used
#'   to calculate confidence intervals for individual studies, see
#'   Details.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.ma The level used to calculate confidence intervals
#'   for meta-analysis estimates.
#' @param fixed A logical indicating whether a fixed effect / common
#'   effect meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects
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
#'   between-study variance tau-squared.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param null.effect A numeric value specifying the effect under the
#'   null hypothesis.
#' @param method.bias A character string indicating which test is to
#'   be used. Either \code{"Begg"}, \code{"Egger"}, or
#'   \code{"Thompson"}, can be abbreviated. See function
#'   \code{\link{metabias}}.
#' @param backtransf A logical indicating whether results for
#'   transformed proportions (argument \code{sm != "PRAW"}) should be
#'   back transformed in printouts and plots. If TRUE (default),
#'   results will be presented as proportions; otherwise transformed
#'   proportions will be shown. See Details for presentation of
#'   confidence intervals.
#' @param pscale A numeric defining a scaling factor for printing of
#'   single event probabilities.
#' @param text.fixed A character string used in printouts and forest
#'   plot to label the pooled fixed effect estimate.
#' @param text.random A character string used in printouts and forest
#'   plot to label the pooled random effects estimate.
#' @param text.predict A character string used in printouts and forest
#'   plot to label the prediction interval.
#' @param text.w.fixed A character string used to label weights of
#'   fixed effect model.
#' @param text.w.random A character string used to label weights of
#'   random effects model.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param subgroup An optional vector to conduct a meta-analysis with
#'   subgroups.
#' @param subgroup.name A character string with a name for the
#'   subgroup variable.
#' @param print.subgroup.name A logical indicating whether the name of
#'   the subgroup variable should be printed in front of the group
#'   labels.
#' @param sep.subgroup A character string defining the separator
#'   between name of subgroup variable and subgroup label.
#' @param test.subgroup A logical value indicating whether to print
#'   results of test for subgroup differences.
#' @param byvar Deprecated argument (replaced by 'subgroup').
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether the addition of
#'   \code{incr} to studies with zero or all events should result in a
#'   warning.
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}} or
#'   \code{\link[metafor]{rma.glmm}}, respectively.
#' @param \dots Additional arguments passed on to
#'   \code{\link[metafor]{rma.glmm}} function.
#' 
#' @details
#' This function provides methods for fixed effect and random effects
#' meta-analysis of single proportions to calculate an overall
#' proportion. Note, you should use R function \code{\link{metabin}}
#' to compare proportions of pairwise comparisons instead of using
#' \code{metaprop} for each treatment arm separately which will break
#' randomisation in randomised controlled trials.
#' 
#' The following transformations of proportions are
#' implemented to calculate an overall proportion:
#' 
#' \itemize{
#' \item Logit transformation (\code{sm = "PLOGIT"}, default)
#' \item Arcsine transformation (\code{sm = "PAS"})
#' \item Freeman-Tukey Double arcsine transformation (\code{sm = "PFT"})
#' \item Log transformation (\code{sm = "PLN"})
#' \item Raw, i.e. untransformed, proportions (\code{sm = "PRAW"})
#' }
#'
#' A generalised linear mixed model (GLMM) - more specific, a random
#' intercept logistic regression model - can be utilised for the
#' meta-analysis of proportions (Stijnen et al., 2010). This is the
#' default method for the logit transformation (argument \code{sm =
#' "PLOGIT"}). Internally, the \code{\link[metafor]{rma.glmm}}
#' function from R package \bold{metafor} is called to fit a GLMM.
#'
#' Classic meta-analysis (Borenstein et al., 2010) utilising the
#' (un)transformed proportions and corresponding standard errors in
#' the inverse variance method is conducted by calling the
#' \code{\link{metagen}} function internally. This is the only
#' available method for all transformations but the logit
#' transformation. The classic meta-analysis model with logit
#' transformed proportions is used by setting argument \code{method =
#' "Inverse"}.
#' 
#' Default settings are utilised for several arguments (assignments
#' using \code{\link{gs}} function). These defaults can be changed for
#' the current R session using the \code{\link{settings.meta}}
#' function.
#' 
#' Furthermore, R function \code{\link{update.meta}} can be used to
#' rerun a meta-analysis with different settings.
#' 
#' \subsection{Choice of transformation / meta-analysis method}{
#' 
#' Contradictory recommendations on the use of transformations of
#' proportions have been published in the literature. For example,
#' Barendregt et al. (2013) recommend the use of the Freeman-Tukey
#' double arcsine transformation instead of the logit transformation
#' whereas Warton & Hui (2011) strongly advise to use generalised
#' linear mixed models with the logit transformation instead of the
#' arcsine transformation.
#'
#' Schwarzer et al. (2019) describe seriously misleading results in a
#' meta-analysis with very different sample sizes due to problems with
#' the back-transformation of the Freeman-Tukey transformation which
#' requires a single sample size (Miller, 1978). Accordingly,
#' Schwarzer et al. (2019) also recommend to use GLMMs for the
#' meta-analysis of single proportions, however, admit that individual
#' study weights are not available with this method. Meta-analysts
#' which require individual study weights should consider the inverse
#' variance method with the arcsine or logit transformation.
#'
#' In order to prevent misleading conclusions for the Freeman-Tukey
#' double arcsine transformation, sensitivity analyses using other
#' transformations or using a range of sample sizes should be
#' conducted (Schwarzer et al., 2019).
#' }
#' 
#' \subsection{Continuity correction}{
#' 
#' If the summary measure is equal to "PLOGIT", "PLN", or "PRAW", a
#' continuity correction is applied if any study has either zero or
#' all events, i.e., an event probability of either 0 or 1.
#'
#' By default, 0.5 is used as continuity correction (argument
#' \code{incr}). This continuity correction is used both to calculate
#' individual study results with confidence limits and to conduct
#' meta-analysis based on the inverse variance method. For GLMMs no
#' continuity correction is used.
#' }
#' 
#' \subsection{Confidence intervals for individual studies}{
#' 
#' Various methods are available to calculate confidence intervals for
#' individual study results (see Agresti & Coull 1998 and Newcombe
#' 1988):
#' \itemize{
#' \item Clopper-Pearson interval also called 'exact' binomial
#'   interval (\code{method.ci = "CP"}, default)
#' \item Wilson Score interval (\code{method.ci = "WS"})
#' \item Wilson Score interval with continuity correction
#'   (\code{method.ci = "WSCC"})
#' \item Agresti-Coull interval (\code{method.ci = "AC"})
#' \item Simple approximation interval (\code{method.ci = "SA"})
#' \item Simple approximation interval with continuity correction
#'   (\code{method.ci = "SACC"})
#' \item Normal approximation interval based on summary measure,
#'   i.e. defined by argument \code{sm} (\code{method.ci = "NAsm"})
#' }
#' 
#' Note, with exception of the normal approximation based on the
#' summary measure, i.e. \code{method.ci = "NAsm"}, the same
#' confidence interval is calculated for individual studies for any
#' summary measure (argument \code{sm}) as only number of events and
#' observations are used in the calculation disregarding the chosen
#' transformation.
#'
#' Results will be presented for transformed proportions if argument
#' \code{backtransf = FALSE}. In this case, argument \code{method.ci =
#' "NAsm"} is used, i.e. confidence intervals based on the normal
#' approximation based on the summary measure.
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
#' from the classic random effects model with the HK method (Knapp and
#' Hartung, 2003; IQWiQ, 2020). An alternative approach is to use the
#' wider confidence interval of classic fixed or random effects
#' meta-analysis and the HK method (Wiksten et al., 2016; Jackson et
#' al., 2017).
#'
#' Argument \code{adhoc.hakn} can be used to choose the \emph{ad hoc}
#' method:
#' \tabular{ll}{
#' \bold{Argument}\tab \bold{\emph{Ad hoc} method} \cr
#' \code{adhoc.hakn = ""}\tab not used \cr
#' \code{adhoc.hakn = "se"}\tab use variance correction if HK standard
#'  error is smaller \cr
#'  \tab than standard error from classic random effects
#'  \cr
#'  \tab meta-analysis (Knapp and Hartung, 2003) \cr
#' \code{adhoc.hakn = "iqwig6"}\tab use variance correction if HK
#'  confidence interval \cr
#'  \tab is narrower than CI from classic random effects model \cr
#'  \tab with DerSimonian-Laird estimator (IQWiG, 2020) \cr
#' \code{adhoc.hakn = "ci"}\tab use wider confidence interval of
#'  classic random effects \cr
#'  \tab and HK meta-analysis \cr
#'  \tab (Hybrid method 2 in Jackson et al., 2017)
#' }
#' 
#' For GLMMs, a method similar to Knapp and Hartung (2003) is
#' implemented, see description of argument \code{tdist} in
#' \code{\link[metafor]{rma.glmm}}, and the \emph{ad hoc} variance
#' correction is not available.
#' }
#' 
#' \subsection{Prediction interval}{
#' 
#' A prediction interval for the proportion in a new study (Higgins et
#' al., 2009) is calculated if arguments \code{prediction} and
#' \code{random} are \code{TRUE}. Note, the definition of
#' prediction intervals varies in the literature. This function
#' implements equation (12) of Higgins et al., (2009) which proposed a
#' \emph{t} distribution with \emph{K-2} degrees of freedom where
#' \emph{K} corresponds to the number of studies in the meta-analysis.
#' }
#'
#' \subsection{Subgroup analysis}{
#' 
#' Argument \code{subgroup} can be used to conduct subgroup analysis for
#' a categorical covariate. The \code{\link{metareg}} function can be
#' used instead for more than one categorical covariate or continuous
#' covariates.
#' }
#' 
#' \subsection{Specify the null hypothesis of test for an overall proportion}{
#'
#' Argument \code{null.effect} can be used to specify the proportion
#' used under the null hypothesis in a test for an overall effect.
#'
#' By default (\code{null.effect = NA}), no hypothesis test is
#' conducted as it is unclear which value is a sensible choice for the
#' data at hand.  An overall proportion of 50\%, for example, could be
#' tested by setting argument \code{null.effect = 0.5}.
#'
#' Note, all tests for an overall effect are two-sided with the
#' alternative hypothesis that the effect is unequal to
#' \code{null.effect}.
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
#' \code{fixed} and \code{random}. Accordingly, the estimate
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"meta"} even if
#' argument \code{random = FALSE}. However, all functions in R
#' package \bold{meta} will adequately consider the values for
#' \code{fixed} and \code{random}. E.g. function
#' \code{\link{print.meta}} will not print results for the random
#' effects model if \code{random = FALSE}.
#' 
#' Argument \code{pscale} can be used to rescale proportions, e.g.
#' \code{pscale = 1000} means that proportions are expressed as events
#' per 1000 observations. This is useful in situations with (very) low
#' event probabilities.
#' }
#' 
#' @return
#' An object of class \code{c("metaprop", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#' \item{event, n, studlab, exclude,}{As defined above.}
#' \item{sm, incr, allincr, addincr, method.ci,}{As defined above.}
#' \item{level, level.ma,}{As defined above.}
#' \item{fixed, random,}{As defined above.}
#' \item{overall, overall.hetstat,}{As defined above.}
#' \item{hakn, adhoc.hakn, method.tau, method.tau.ci,}{As defined above.}
#' \item{tau.preset, TE.tau, null.hypothesis,}{As defined above.}
#' \item{method.bias, tau.common, title, complab, outclab,}{As defined
#'   above.}
#' \item{subgroup, subgroup.name, print.subgroup.name, sep.subgroup, warn}{As defined
#'   above.}
#' \item{TE, seTE}{Estimated (un)transformed proportion and its
#'   standard error for individual studies.}
#' \item{lower, upper}{Lower and upper confidence interval limits for
#'   individual studies.}
#' \item{zval, pval}{z-value and p-value for test of treatment effect
#'   for individual studies.}
#' \item{w.fixed, w.random}{Weight of individual studies (in fixed and
#'   random effects model).}
#' \item{TE.fixed, seTE.fixed}{Estimated overall (un)transformed
#'   proportion and standard error (fixed effect model).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (fixed effect model).}
#' \item{statistic.fixed, pval.fixed}{z-value and p-value for test of
#'   overall effect (fixed effect model).}
#' \item{TE.random, seTE.random}{Estimated overall (un)transformed
#'   proportion and standard error (random effects model).}
#' \item{lower.random, upper.random}{Lower and upper confidence
#'   interval limits (random effects model).}
#' \item{statistic.random, pval.random}{z-value or t-value and
#'   corresponding p-value for test of overall effect (random effects
#'   model).}
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
#' \item{method}{A character string indicating method used for
#'   pooling: \code{"Inverse"}}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn=TRUE}).}
#' \item{bylevs}{Levels of grouping variable - if \code{subgroup} is not
#'   missing.}
#' \item{TE.fixed.w, seTE.fixed.w}{Estimated treatment effect and
#'   standard error in subgroups (fixed effect model) - if
#'   \code{subgroup} is not missing.}
#' \item{lower.fixed.w, upper.fixed.w}{Lower and upper confidence
#'   interval limits in subgroups (fixed effect model) - if
#'   \code{subgroup} is not missing.}
#' \item{statistic.fixed.w, pval.fixed.w}{z-value and p-value for test
#'   of treatment effect in subgroups (fixed effect model) - if
#'   \code{subgroup} is not missing.}
#' \item{TE.random.w, seTE.random.w}{Estimated treatment effect and
#'   standard error in subgroups (random effects model) - if
#'   \code{subgroup} is not missing.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{subgroup} is not missing.}
#' \item{statistic.random.w, pval.random.w}{z-value or t-value and
#'   corresponding p-value for test of treatment effect in subgroups
#'   (random effects model) - if \code{subgroup} is not missing.}
#' \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
#'   random effects model) - if \code{subgroup} is not missing.}
#' \item{df.hakn.w}{Degrees of freedom for test of treatment effect
#'   for Hartung-Knapp method in subgroups - if \code{subgroup} is not
#'   missing and \code{hakn=TRUE}.}
#' \item{n.harmonic.mean.w}{Harmonic mean of number of observations in
#'   subgroups (for back transformation of Freeman-Tukey Double
#'   arcsine transformation) - if \code{subgroup} is not missing.}
#' \item{event.w}{Number of events in subgroups - if \code{subgroup} is
#'   not missing.}
#' \item{n.w}{Number of observations in subgroups - if \code{subgroup} is
#'   not missing.}
#' \item{k.w}{Number of studies combined within subgroups - if
#'   \code{subgroup} is not missing.}
#' \item{k.all.w}{Number of all studies in subgroups - if \code{subgroup}
#'   is not missing.}
#' \item{Q.w.fixed}{Overall within subgroups heterogeneity statistic Q
#'   (based on fixed effect model) - if \code{subgroup} is not missing.}
#' \item{Q.w.random}{Overall within subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{subgroup} is not
#'   missing (only calculated if argument \code{tau.common} is TRUE).}
#' \item{df.Q.w}{Degrees of freedom for test of overall within
#'   subgroups heterogeneity - if \code{subgroup} is not missing.}
#' \item{pval.Q.w.fixed}{P-value of within subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{subgroup} is
#'   not missing.}
#' \item{pval.Q.w.random}{P-value of within subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{subgroup} is
#'   not missing.}
#' \item{Q.b.fixed}{Overall between subgroups heterogeneity statistic
#'   Q (based on fixed effect model) - if \code{subgroup} is not
#'   missing.}
#' \item{Q.b.random}{Overall between subgroups heterogeneity statistic
#'   Q (based on random effects model) - if \code{subgroup} is not
#'   missing.}
#' \item{df.Q.b}{Degrees of freedom for test of overall between
#'   subgroups heterogeneity - if \code{subgroup} is not missing.}
#' \item{pval.Q.b.fixed}{P-value of between subgroups heterogeneity
#'   statistic Q (based on fixed effect model) - if \code{subgroup} is
#'   not missing.}
#' \item{pval.Q.b.random}{P-value of between subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{subgroup} is
#'   not missing.}
#' \item{tau.w}{Square-root of between-study variance within subgroups
#'   - if \code{subgroup} is not missing.}
#' \item{H.w}{Heterogeneity statistic H within subgroups - if
#'   \code{subgroup} is not missing.}
#' \item{lower.H.w, upper.H.w}{Lower and upper confidence limit for
#'   heterogeneity statistic H within subgroups - if \code{subgroup} is
#'   not missing.}
#' \item{I2.w}{Heterogeneity statistic I\eqn{^2} within subgroups - if
#'   \code{subgroup} is not missing.}
#' \item{lower.I2.w, upper.I2.w}{Lower and upper confidence limit for
#'   heterogeneity statistic I\eqn{^2} within subgroups - if \code{subgroup} is
#'   not missing.}
#' \item{incr.event}{Increment added to number of events.}
#' \item{keepdata}{As defined above.}
#' \item{data}{Original data (set) used in function call (if
#'   \code{keepdata=TRUE}).}
#' \item{subset}{Information on subset of original data used in
#'   meta-analysis (if \code{keepdata=TRUE}).}
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
#' @seealso \code{\link{update.meta}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{print.meta}},
#'   \code{\link{forest.meta}}
#' 
#' @references
#' Agresti A & Coull BA (1998):
#' Approximate is better than "exact" for interval estimation of
#' binomial proportions.
#' \emph{The American Statistician},
#' \bold{52}, 119--26
#' 
#' Barendregt JJ, Doi SA, Lee YY, Norman RE, Vos T (2013):
#' Meta-analysis of prevalence.
#' \emph{Journal of Epidemiology and Community Health},
#' \bold{67}, 974--8
#' 
#' Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2010):
#' A basic introduction to fixed-effect and random-effects models for
#' meta-analysis.
#' \emph{Research Synthesis Methods},
#' \bold{1}, 97--111
#' 
#' DerSimonian R & Laird N (1986):
#' Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials},
#' \bold{7}, 177--88
#' 
#' Edward JM et al. (2006):
#' Adherence to antiretroviral therapy in sub-saharan
#' Africa and North America - a meta-analysis.
#' \emph{Journal of the American Medical Association},
#' \bold{296}, 679--90
#' 
#' Freeman MF & Tukey JW (1950):
#' Transformations related to the angular and the square root.
#' \emph{Annals of Mathematical Statistics},
#' \bold{21}, 607--11
#' 
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009):
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{172}, 137--59
#' 
#' Hartung J, Knapp G (2001a):
#' On tests of the overall treatment effect in meta-analysis with
#' normally distributed responses.
#' \emph{Statistics in Medicine},
#' \bold{20}, 1771--82
#' 
#' Hartung J, Knapp G (2001b):
#' A refined method for the meta-analysis of controlled clinical
#' trials with binary outcome.
#' \emph{Statistics in Medicine},
#' \bold{20}, 3875--89
#'
#' IntHout J, Ioannidis JPA, Borm GF (2014):
#' The Hartung-Knapp-Sidik-Jonkman method for random effects
#' meta-analysis is straightforward and considerably outperforms the
#' standard DerSimonian-Laird method.
#' \emph{BMC Medical Research Methodology},
#' \bold{14}, 25
#' 
#' IQWiG (2020):
#' General Methods: Version 6.0.
#' \url{https://www.iqwig.de/en/about-us/methods/methods-paper/}
#'
#' Jackson D, Law M, Rücker G, Schwarzer G (2017): 
#' The Hartung-Knapp modification for random-effects meta-analysis: A
#' useful refinement but are there any residual concerns?
#' \emph{Statistics in Medicine},
#' \bold{36}, 3923--34
#' 
#' Langan D, Higgins JPT, Jackson D, Bowden J, Veroniki AA,
#' Kontopantelis E, et al. (2019):
#' A comparison of heterogeneity variance estimators in simulated
#' random-effects meta-analyses.
#' \emph{Research Synthesis Methods},
#' \bold{10}, 83--98
#' 
#' Miller JJ (1978):
#' The inverse of the Freeman-Tukey double arcsine transformation.
#' \emph{The American Statistician},
#' \bold{32}, 138
#' 
#' Newcombe RG (1998):
#' Two-sided confidence intervals for the single proportion:
#' comparison of seven methods.
#' \emph{Statistics in Medicine},
#' \bold{17}, 857--72
#' 
#' Pettigrew HM, Gart JJ, Thomas DG (1986):
#' The bias and higher cumulants of the logarithm of a binomial
#' variate.
#' \emph{Biometrika},
#' \bold{73}, 425--35
#' 
#' Schwarzer G, Chemaitelly H, Abu-Raddad LJ, Rücker G (2019):
#' Seriously misleading results using inverse of Freeman-Tukey double
#' arcsine transformation in meta-analysis of single proportions.
#' \emph{Research Synthesis Methods},
#' \bold{10}, 476--83
#' 
#' Stijnen T, Hamza TH, Ozdemir P (2010):
#' Random effects meta-analysis of event outcome in the framework of
#' the generalized linear mixed model with applications in sparse
#' data.
#' \emph{Statistics in Medicine},
#' \bold{29}, 3046--67
#' 
#' Viechtbauer W (2010):
#' Conducting meta-analyses in R with the metafor package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#'
#' Warton DI, Hui FKC (2011):
#' The arcsine is asinine: the analysis of proportions in ecology.
#' \emph{Ecology},
#' \bold{92}, 3--10
#'
#' Wiksten A, Rücker G, Schwarzer G (2016):
#' Hartung-Knapp method is not always conservative compared with
#' fixed-effect meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{35}, 2503--15
#' 
#' @examples
#' # Meta-analysis using generalised linear mixed model
#' #
#' metaprop(4:1, 10 * 1:4)
#' 
#' # Apply various classic meta-analysis methods to estimate
#' # proportions
#' #
#' m1 <- metaprop(4:1, 10 * 1:4, method = "Inverse")
#' m2 <- update(m1, sm = "PAS")
#' m3 <- update(m1, sm = "PRAW")
#' m4 <- update(m1, sm = "PLN")
#' m5 <- update(m1, sm = "PFT")
#' #
#' m1
#' m2
#' m3
#' m4
#' m5
#' #
#' forest(m1)
#' \dontrun{
#' forest(m2)
#' forest(m3)
#' forest(m3, pscale = 100)
#' forest(m4)
#' forest(m5)
#' }
#' 
#' # Do not back transform results, e.g. print logit transformed
#' # proportions if sm = "PLOGIT" and store old settings
#' #
#' oldset <- settings.meta(backtransf = FALSE)
#' #
#' m6  <- metaprop(4:1, c(10, 20, 30, 40), method = "Inverse")
#' m7  <- update(m6, sm = "PAS")
#' m8  <- update(m6, sm = "PRAW")
#' m9  <- update(m6, sm = "PLN")
#' m10 <- update(m6, sm = "PFT")
#' #
#' forest(m6)
#' \dontrun{
#' forest(m7)
#' forest(m8)
#' forest(m8, pscale = 100)
#' forest(m9)
#' forest(m10)
#' }
#' 
#' # Use old settings
#' #
#' settings.meta(oldset)
#' 
#' # Examples with zero events
#' #
#' m1 <- metaprop(c(0, 0, 10, 10), rep(100, 4), method = "Inverse")
#' m2 <- metaprop(c(0, 0, 10, 10), rep(100, 4), incr = 0.1, method = "Inverse")
#' #
#' m1
#' m2
#' #
#' \dontrun{
#' forest(m1)
#' forest(m2)
#' }
#' 
#' # Example from Miller (1978):
#' #
#' death <- c(3, 6, 10, 1)
#' animals <- c(11, 17, 21, 6)
#' #
#' m3 <- metaprop(death, animals, sm = "PFT")
#' forest(m3)
#' 
#' # Data examples from Newcombe (1998)
#' # - apply various methods to estimate confidence intervals for
#' #   individual studies
#' #
#' event <- c(81, 15, 0, 1)
#' n <- c(263, 148, 20, 29)
#' #
#' m1 <- metaprop(event, n, method.ci = "SA", method = "Inverse")
#' m2 <- update(m1, method.ci = "SACC")
#' m3 <- update(m1, method.ci = "WS")
#' m4 <- update(m1, method.ci = "WSCC")
#' m5 <- update(m1, method.ci = "CP")
#' #
#' lower <- round(rbind(NA, m1$lower, m2$lower, NA, m3$lower,
#'                      m4$lower, NA, m5$lower), 4)
#' upper <- round(rbind(NA, m1$upper, m2$upper, NA, m3$upper,
#'                      m4$upper, NA, m5$upper), 4)
#' #
#' tab1 <- data.frame(
#'   scen1 = meta:::formatCI(lower[, 1], upper[, 1]),
#'   scen2 = meta:::formatCI(lower[, 2], upper[, 2]),
#'   scen3 = meta:::formatCI(lower[, 3], upper[, 3]),
#'   scen4 = meta:::formatCI(lower[, 4], upper[, 4]),
#'   stringsAsFactors = FALSE
#'   )
#' names(tab1) <- c("r=81, n=263", "r=15, n=148",
#'                  "r=0, n=20", "r=1, n=29")
#' row.names(tab1) <- c("Simple", "- SA", "- SACC",
#'                      "Score", "- WS", "- WSCC",
#'                      "Binomial", "- CP")
#' tab1[is.na(tab1)] <- ""
#' # Newcombe (1998), Table I, methods 1-5:
#' tab1
#' 
#' # Same confidence interval, i.e. unaffected by choice of summary
#' # measure
#' #
#' print(metaprop(event, n, method.ci = "WS", method = "Inverse"), ma = FALSE)
#' print(metaprop(event, n, sm = "PLN", method.ci = "WS"), ma = FALSE)
#' print(metaprop(event, n, sm = "PFT", method.ci = "WS"), ma = FALSE)
#' print(metaprop(event, n, sm = "PAS", method.ci = "WS"), ma = FALSE)
#' print(metaprop(event, n, sm = "PRAW", method.ci = "WS"), ma = FALSE)
#' 
#' # Different confidence intervals as argument sm = "NAsm"
#' #
#' print(metaprop(event, n, method.ci = "NAsm", method = "Inverse"), ma = FALSE)
#' print(metaprop(event, n, sm = "PLN", method.ci = "NAsm"), ma = FALSE)
#' print(metaprop(event, n, sm = "PFT", method.ci = "NAsm"), ma = FALSE)
#' print(metaprop(event, n, sm = "PAS", method.ci = "NAsm"), ma = FALSE)
#' print(metaprop(event, n, sm = "PRAW", method.ci = "NAsm"), ma = FALSE)
#' 
#' # Different confidence intervals as argument backtransf = FALSE.
#' # Accordingly, method.ci = "NAsm" used internally.
#' #
#' print(metaprop(event, n, method.ci = "WS", method = "Inverse"),
#'       ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PLN", method.ci = "WS"),
#'       ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PFT", method.ci = "WS"),
#'       ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PAS", method.ci = "WS"),
#'       ma = FALSE, backtransf = FALSE)
#' print(metaprop(event, n, sm = "PRAW", method.ci = "WS"),
#'       ma = FALSE, backtransf = FALSE)
#' 
#' # Same results (printed on original and log scale, respectively)
#' #
#' print(metaprop(event, n, sm = "PLN", method.ci = "NAsm"), ma = FALSE)
#' print(metaprop(event, n, sm = "PLN"), ma = FALSE, backtransf = FALSE)
#' # Results for first study (on log scale)
#' round(log(c(0.3079848, 0.2569522, 0.3691529)), 4)
#' 
#' # Print results as events per 1000 observations
#' #
#' print(metaprop(6:8, c(100, 1200, 1000), method = "Inverse"),
#'       pscale = 1000, digits = 1)
#' 
#' @export metaprop


metaprop <- function(event, n, studlab,
                     ##
                     data = NULL, subset = NULL, exclude = NULL,
                     method,
                     ##
                     sm = gs("smprop"),
                     ##
                     incr = gs("incr"), allincr = gs("allincr"),
                     addincr = gs("addincr"),
                     ##
                     method.ci = gs("method.ci.prop"),
                     level = gs("level"), level.ma = gs("level.ma"),
                     fixed = gs("fixed"),
                     random = gs("random") | !is.null(tau.preset),
                     overall = fixed | random,
                     overall.hetstat = fixed | random,
                     ##
                     hakn = gs("hakn"), adhoc.hakn = gs("adhoc.hakn"),
                     method.tau,
                     method.tau.ci = gs("method.tau.ci"),
                     tau.preset = NULL, TE.tau = NULL,
                     tau.common = gs("tau.common"),
                     ##
                     prediction = gs("prediction"),
                     level.predict = gs("level.predict"),
                     ##
                     null.effect = NA,
                     ##
                     method.bias = gs("method.bias"),
                     ##
                     backtransf = gs("backtransf"),
                     pscale = 1,
                     ##
                     text.fixed = gs("text.fixed"),
                     text.random = gs("text.random"),
                     text.predict = gs("text.predict"),
                     text.w.fixed = gs("text.w.fixed"),
                     text.w.random = gs("text.w.random"),
                     ##
                     title = gs("title"), complab = gs("complab"),
                     outclab = "",
                     ##
                     subgroup, subgroup.name = NULL,
                     print.subgroup.name = gs("print.subgroup.name"),
                     sep.subgroup = gs("sep.subgroup"),
                     test.subgroup = gs("test.subgroup"),
                     byvar,
                     ##
                     keepdata = gs("keepdata"),
                     warn = gs("warn"), warn.deprecated = gs("warn.deprecated"),
                     ##
                     control = NULL,
                     ...
                     ) {
  
  
  ##
  ##
  ## (1) Check and set arguments
  ##
  ##
  if (missing(method))
    method <- if (sm == "PLOGIT") "GLMM" else "Inverse"
  else
    method <- setchar(method, .settings$meth4prop)
  is.glmm <- method == "GLMM"
  ##
  chknull(sm)
  sm <- setchar(sm, .settings$sm4prop)
  ##
  chklevel(level)
  ##
  chklogical(hakn)
  missing.adhoc.hakn <- missing(adhoc.hakn)
  adhoc.hakn <- setchar(adhoc.hakn, .settings$adhoc4hakn)
  if (missing(method.tau))
    method.tau <- if (method == "GLMM") "ML" else gs("method.tau")
  method.tau <- setchar(method.tau, .settings$meth4tau)
  if (is.null(method.tau.ci))
    method.tau.ci <- if (method.tau == "DL") "J" else "QP"
  method.tau.ci <- setchar(method.tau.ci, .settings$meth4tau.ci)
  chklogical(tau.common)
  ##
  chklogical(prediction)
  chklevel(level.predict)
  ##
  if (!anyNA(null.effect) | length(null.effect) != 1)
    chknumeric(null.effect, min = 0, max = 1, length = 1)
  ##
  method.bias <- setmethodbias(method.bias)
  ##
  chklogical(backtransf)
  ##
  chknumeric(pscale, length = 1)
  if (!backtransf & pscale != 1 & !is.untransformed(sm)) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  ##
  if (!is.null(text.fixed))
    chkchar(text.fixed, length = 1)
  if (!is.null(text.random))
    chkchar(text.random, length = 1)
  if (!is.null(text.predict))
    chkchar(text.predict, length = 1)
  if (!is.null(text.w.fixed))
    chkchar(text.w.fixed, length = 1)
  if (!is.null(text.w.random))
    chkchar(text.w.random, length = 1)
  ##
  chklogical(keepdata)
  ##
  ## Additional arguments / checks for metainc objects
  ##
  fun <- "metaprop"
  ##
  chklogical(allincr)
  chklogical(addincr)
  method.ci <- setchar(method.ci, .settings$ci4prop)
  chklogical(warn)
  ##
  if (is.glmm & sm != "PLOGIT")
    stop("Generalised linear mixed models only possible with ",
         "argument 'sm = \"PLOGIT\"'.")
  ##
  if (is.glmm & method.tau != "ML")
    stop("Generalised linear mixed models only possible with ",
         "argument 'method.tau = \"ML\"'.")
   ##
  if (is.glmm & hakn & adhoc.hakn != "") {
    if (!missing.adhoc.hakn)
      warning("Ad hoc variance correction for Hartung-Knapp method ",
              "not available for GLMMs.",
              call. = FALSE)
    adhoc.hakn <- ""
  }
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  chklogical(warn.deprecated)
  ##
  level.ma <- deprecated(level.ma, missing(level.ma), args, "level.comb",
                         warn.deprecated)
  chklevel(level.ma)
  ##
  fixed <- deprecated(fixed, missing(fixed), args, "comb.fixed",
                      warn.deprecated)
  chklogical(fixed)
  ##
  random <- deprecated(random, missing(random), args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  ##
  missing.subgroup.name <- missing(subgroup.name)
  subgroup.name <-
    deprecated(subgroup.name, missing.subgroup.name, args, "bylab",
               warn.deprecated)
  ##
  print.subgroup.name <-
    deprecated(print.subgroup.name, missing(print.subgroup.name),
               args, "print.byvar", warn.deprecated)
  print.subgroup.name <- replaceNULL(print.subgroup.name, FALSE)
  chklogical(print.subgroup.name)
  ##
  sep.subgroup <-
    deprecated(sep.subgroup, missing(sep.subgroup), args, "byseparator",
               warn.deprecated)
  if (!is.null(sep.subgroup))
    chkchar(sep.subgroup, length = 1)
  ##
  ## Some more checks
  ##
  chklogical(overall)
  chklogical(overall.hetstat)
  
  
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
  ## Catch 'event' and 'n' from data:
  ##
  event <- eval(mf[[match("event", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  chknull(event)
  k.All <- length(event)
  ##
  n <- eval(mf[[match("n", names(mf))]],
            data, enclos = sys.frame(sys.parent()))
  chknull(n)
  ##
  ## Catch 'incr' from data:
  ##
  if (!missing(incr))
    incr <- eval(mf[[match("incr", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  chknumeric(incr, min = 0)
  ##
  ## Catch 'studlab', 'subgroup', 'subset', and 'exclude' from data:
  ##
  studlab <- eval(mf[[match("studlab", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  studlab <- setstudlab(studlab, k.All)
  ##
  missing.subgroup <- missing(subgroup)
  subgroup <- eval(mf[[match("subgroup", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
  missing.byvar <- missing(byvar)
  byvar <- eval(mf[[match("byvar", names(mf))]],
                data, enclos = sys.frame(sys.parent()))
  ##
  subgroup <- deprecated2(subgroup, missing.subgroup, byvar, missing.byvar,
                          warn.deprecated)
  by <- !is.null(subgroup)
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
  chklength(n, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (length(incr) > 1)
    chklength(incr, k.All, fun)
  ##
  if (by) {
    chklength(subgroup, k.All, fun)
    chklogical(test.subgroup)
  }
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
    warning("Value for argument 'tau.common' set to FALSE as argument 'subgroup' is missing.")
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
  ##
  ## (5) Store complete dataset in list object data
  ##     (if argument keepdata is TRUE)
  ##
  ##
  if (keepdata) {
    if (nulldata)
      data <- data.frame(.event = event)
    else
      data$.event <- event
    ##
    data$.n <- n
    data$.studlab <- studlab
    ##
    data$.incr <- incr
    ##
    if (by)
      data$.subgroup <- subgroup
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
  }
  
  
  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    event <- event[subset]
    n <- n[subset]
    studlab <- studlab[subset]
    ##
    exclude <- exclude[subset]
    ##
    if (length(incr) > 1)
      incr <- incr[subset]
    ##
    if (by)
      subgroup <- subgroup[subset]
  }
  ##
  ## Determine total number of studies
  ##
  k.all <- length(event)
  ##
  if (k.all == 0)
    stop("No studies to combine in meta-analysis.")
  ##
  ## No meta-analysis for a single study
  ##
  if (k.all == 1) {
    fixed  <- FALSE
    random <- FALSE
    prediction  <- FALSE
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  ##
  ## Check variable values
  ##
  chknumeric(event, 0)
  chknumeric(n, 0, zero = TRUE)
  ##
  if (any(event > n, na.rm = TRUE))
    stop("Number of events must not be larger than number of observations")
  ##
  ## Recode integer as numeric:
  ##
  event <- int2num(event)
  n     <- int2num(n)
  ##
  ## Check for whole numbers
  ##
  if (method.ci != "NAsm") {
    if (any(!is.wholenumber(event), na.rm = TRUE)) {
      warning("Normal approximation confidence interval ",
              "(argument method.ci = \"NAsm\") used as\n",
              "at least one number of events contains a non-integer value.",
              call. = FALSE)
      method.ci <- "NAsm"
    }
    else if (any(!is.wholenumber(n), na.rm = TRUE)) {
      warning("Normal approximation confidence interval ",
              "(argument method.ci = \"NAsm\") used as\n",
              "at least one sample size contains a non-integer value.",
              call. = FALSE)
      method.ci <- "NAsm"
    }
  }
  ##
  if (by) {
    chkmiss(subgroup)
    ##
    if (missing.subgroup.name & is.null(subgroup.name)) {
      if (!missing.subgroup)
        subgroup.name <- byvarname(mf[[match("subgroup", names(mf))]])
      else if (!missing.byvar)
        subgroup.name <- byvarname(mf[[match("byvar", names(mf))]])
    }
  }
  ##
  if (!is.null(subgroup.name))
    chkchar(subgroup.name, length = 1)
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  sel <- switch(sm,
                PLOGIT = event == 0 | (n - event) == 0,
                PAS = rep(FALSE, length(event)),
                PFT = rep(FALSE, length(event)),
                PLN =    event == 0 | (n - event) == 0,
                PRAW =   event == 0 | (n - event) == 0)
  ##
  sparse <- any(sel, na.rm = TRUE)
  ##
  ## No need to add anything to cell counts for arcsine transformation
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
  if (sm == "PLOGIT") {
    TE <- log((event + incr.event) / (n - event + incr.event))
    seTE <- sqrt(1 / (event + incr.event) +
                 1 / ((n - event + incr.event)))
    transf.null.effect <- log(null.effect / (1 - null.effect))
  }
  else if (sm == "PAS") {
    TE <- asin(sqrt(event / n))
    seTE <- sqrt(1 / (4 * n))
    transf.null.effect <- asin(sqrt(null.effect))
  }
  else if (sm == "PFT") {
    TE <- 0.5 * (asin(sqrt(event / (n + 1))) + asin(sqrt((event + 1) / (n + 1))))
    seTE <- sqrt(1 / (4 * n + 2))
    transf.null.effect <- asin(sqrt(null.effect))
  }
  else if (sm == "PLN") {
    TE <- log((event + incr.event) / (n + incr.event))
    ## Hartung, Knapp (2001), p. 3880, formula (18):
    seTE <- ifelse(event == n,
                   sqrt(1 / event                - 1 / (n + incr.event)),
                   sqrt(1 / (event + incr.event) - 1 / (n + incr.event))
                   )
    transf.null.effect <- log(null.effect)
  }
  else if (sm == "PRAW") {
    TE <- event / n
    seTE <- sqrt((event + incr.event) * (n - event + incr.event) /
                 (n + 2 * incr.event)^3)
    transf.null.effect <- null.effect
  }
  ##
  ## Calculate confidence intervals
  ##
  NAs <- rep(NA, k.all)
  ##
  if (method.ci == "CP") {
    lower.study <- upper.study <- NAs
    for (i in 1:k.all) {
      if (!is.na(event[i] & !is.na(n[i]))) {
        cint <- binom.test(event[i], n[i], conf.level = level)
        ##
        lower.study[i] <- cint$conf.int[[1]]
        upper.study[i] <- cint$conf.int[[2]]
      }
      else {
        lower.study[i] <- NA
        upper.study[i] <- NA
      }
    }
  }
  ##
  else {
    if (method.ci == "WS")
      ci.study <- ciWilsonScore(event, n, level = level)
    ##
    else if (method.ci == "WSCC")
      ci.study <- ciWilsonScore(event, n, level = level, correct = TRUE)
    ##
    else if (method.ci == "AC")
      ci.study <- ciAgrestiCoull(event, n, level = level)
    ##
    else if (method.ci == "SA")
      ci.study <- ciSimpleAsymptotic(event, n, level = level)
    ##
    else if (method.ci == "SACC")
      ci.study <- ciSimpleAsymptotic(event, n, level = level, correct = TRUE)
    ##
    else if (method.ci == "NAsm")
      ci.study <- ci(TE, seTE, level = level)
    ##
    lower.study <- ci.study$lower
    upper.study <- ci.study$upper
  }
  ##  
  if (method.ci == "NAsm") {
    if (sm == "PLOGIT") {
      lower.study <- logit2p(lower.study)
      upper.study <- logit2p(upper.study)
    }
    ##
    else if (sm == "PAS") {
      lower.study <- asin2p(lower.study, value = "lower", warn = FALSE)
      upper.study <- asin2p(upper.study, value = "upper", warn = FALSE)
    }
    ##
    else if (sm == "PFT") {
      lower.study <- asin2p(lower.study, n, value = "lower", warn = FALSE)
      upper.study <- asin2p(upper.study, n, value = "upper", warn = FALSE)
    }
    ##
    else if (sm == "PLN") {
      lower.study <- exp(lower.study)
      upper.study <- exp(upper.study)
    }
    ##
    lower.study[lower.study < 0] <- 0
    upper.study[upper.study > 1] <- 1
  }
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(event[!exclude]) & !is.na(n[!exclude]))
  ##
  if (k == 1 & hakn)
    hakn <- FALSE
  ##
  if (is.glmm & k > 0) {
    glmm.fixed <- rma.glmm(xi = event[!exclude], ni = n[!exclude],
                           method = "FE", test = ifelse(hakn, "t", "z"),
                           level = 100 * level.ma,
                           measure = "PLO", control = control,
                           ...)
    ##
    TE.fixed   <- as.numeric(glmm.fixed$b)
    seTE.fixed <- as.numeric(glmm.fixed$se)
    ##
    w.fixed <- rep(NA, length(event))
  }
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
               ##
               sm = sm,
               level = level,
               level.ma = level.ma,
               fixed = fixed,
               random = random,
               overall = overall,
               overall.hetstat = overall.hetstat,
               ##
               hakn = hakn, adhoc.hakn = adhoc.hakn,
               method.tau = method.tau, method.tau.ci = method.tau.ci,
               tau.preset = tau.preset,
               TE.tau = TE.tau,
               tau.common = FALSE,
               ##
               prediction = prediction,
               level.predict = level.predict,
               ##
               null.effect = transf.null.effect,
               ##
               method.bias = method.bias,
               ##
               backtransf = backtransf,
               ##
               text.fixed = text.fixed, text.random = text.random,
               text.predict = text.predict,
               text.w.fixed = text.w.fixed, text.w.random = text.w.random,
               ##
               title = title, complab = complab, outclab = outclab,
               ##
               keepdata = FALSE,
               warn = warn,
               ##
               control = control)
  ##
  if (by & tau.common & !is.glmm) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau, "",
                   TE.tau, level.ma, subgroup, control)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(event = event, n = n,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              sparse = sparse,
              allincr = allincr, addincr = addincr,
              method.ci = method.ci,
              incr.event = incr.event)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$n.e <- NULL
  m$n.c <- NULL
  m$label.e <- ""
  m$label.c <- ""
  m$label.left <- ""
  m$label.right <- ""
  ##
  res <- c(res, m)
  res$null.effect <- null.effect
  ##
  ## Add data
  ##
  if (is.glmm & k > 0) {
    ##
    ci.f <- ci(TE.fixed, seTE.fixed, level = level.ma,
               null.effect = transf.null.effect)
    ##
    res$method <- "GLMM"
    ##
    res$TE.fixed <- TE.fixed
    res$seTE.fixed <- seTE.fixed
    res$w.fixed <- w.fixed
    res$lower.fixed <- ci.f$lower
    res$upper.fixed <- ci.f$upper
    res$statistic.fixed <- ci.f$statistic
    res$pval.fixed <- ci.f$p
    res$zval.fixed <- ci.f$statistic
    ##
    if (sum(!exclude) > 1 &
        sum(event[!exclude], na.rm = TRUE) > 0 &
        any(event[!exclude] != n[!exclude]))
      glmm.random <- rma.glmm(xi = event[!exclude], ni = n[!exclude],
                              method = method.tau,
                              test = ifelse(hakn, "t", "z"),
                              level = 100 * level.ma,
                                measure = "PLO", control = control,
                              ...)
    else {
      ##
      ## Fallback to fixed effect model due to small number of studies
      ## or zero or all events
      ##
      glmm.random <- glmm.fixed
    }
    ##
    TE.random   <- as.numeric(glmm.random$b)
    seTE.random <- as.numeric(glmm.random$se)
    ##
    ci.r <- ci(TE.random, seTE.random, level = level.ma,
               null.effect = transf.null.effect, df = if (hakn) k - 1)
    ##
    res$w.random <- rep(NA, length(event))
    ##
    res$TE.random <- TE.random
    res$seTE.random <- seTE.random
    res$lower.random <- ci.r$lower
    res$upper.random <- ci.r$upper
    res$statistic.random <- ci.r$statistic
    res$pval.random <- ci.r$p
    res$zval.random <- ci.r$statistic
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
    H <- calcH(res$Q, res$df.Q, level.ma)
    res$H <- H$TE
    res$lower.H <- H$lower
    res$upper.H <- H$upper
    ##
    I2 <- isquared(res$Q, res$df.Q, level.ma)
    res$I2 <- I2$TE
    res$lower.I2 <- I2$lower
    res$upper.I2 <- I2$upper
    ##
    res$Rb <- NA
    res$lower.Rb <- NA
    res$upper.Rb <- NA
    ##
    res$.glmm.fixed  <- glmm.fixed
    res$.glmm.random <- glmm.random
    res$version.metafor <- packageDescription("metafor")$Version
    ##
    if (by) {
      n.by <- length(unique(subgroup[!exclude]))
      if (n.by > 1)
        subgroup.glmm <- factor(subgroup[!exclude], bylevs(subgroup[!exclude]))
      ##
      glmm.random.by <-
        try(suppressWarnings(rma.glmm(xi = event[!exclude], ni = n[!exclude],
                                      mods =
                                        if (n.by > 1) ~ subgroup.glmm else NULL,
                                      method = method.tau,
                                      test = ifelse(hakn, "t", "z"),
                                      level = 100 * level.ma,
                                      measure = "PLO", control = control,
                                      ...)),
            silent = TRUE)
      ##
      if ("try-error" %in% class(glmm.random.by))
        if (grepl(paste0("Number of parameters to be estimated is ",
                         "larger than the number of observations"),
                  glmm.random.by)) {
          glmm.random.by <-
            suppressWarnings(rma.glmm(xi = event[!exclude], ni = n[!exclude],
                                      mods =
                                        if (n.by > 1) ~ subgroup.glmm else NULL,
                                      method = "FE",
                                      test = ifelse(hakn, "t", "z"),
                                      level = 100 * level.ma,
                                      measure = "PLO", control = control,
                                      ...))
        }
        else
          stop(glmm.random.by)
      ##
      Q.r <- glmm.random.by$QE.Wld
      df.Q.r <- glmm.random.by$k - glmm.random.by$p
      ##
      H.r  <- calcH(Q.r, df.Q.r, level.ma)
      I2.r <- isquared(Q.r, df.Q.r, level.ma)
      ##
      hcc <- list(tau2.resid = glmm.random.by$tau2,
                  lower.tau2.resid = NA,
                  upper.tau2.resid = NA,
                  ##
                  tau.resid = sqrt(glmm.random.by$tau2),
                  lower.tau.resid = NA,
                  upper.tau.resid = NA,
                  sign.lower.tau.resid = "",
                  sign.upper.tau.resid = "",
                  ##
                  Q.resid = Q.r,
                  df.Q.resid = df.Q.r,
                  pval.Q.resid = pvalQ(Q.r, df.Q.r),
                  ##
                  H.resid = H.r$TE,
                  lower.H.resid = H.r$lower,
                  upper.H.resid = H.r$upper,
                  ##
                  I2.resid = I2.r$TE,
                  lower.I2.resid = I2.r$lower,
                  upper.I2.resid = I2.r$upper
                  )
    }
  }
  ##
  res$lower <- lower.study
  res$upper <- upper.study
  ##
  res$pscale <- pscale
  ##
  res$call <- match.call()
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
    res$subgroup <- subgroup
    res$subgroup.name <- subgroup.name
    res$print.subgroup.name <- print.subgroup.name
    res$sep.subgroup <- sep.subgroup
    res$test.subgroup <- test.subgroup
    res$tau.common <- tau.common
    ##
    if (!tau.common)
      res <- c(res, subgroup(res))
    else if (!is.null(tau.preset))
      res <- c(res, subgroup(res, tau.preset))
    else {
      if (is.glmm)
        res <- c(res, subgroup(res, NULL,
                               factor(res$subgroup, bylevs(res$subgroup)), ...))
      else
        res <- c(res, subgroup(res, hcc$tau.resid))
    }
    ##
    if (!tau.common || !is.null(tau.preset)) {
      res$tau2.resid <- res$lower.tau2.resid <- res$upper.tau2.resid <- NA
      res$tau.resid <- res$lower.tau.resid <- res$upper.tau.resid <- NA
      ##
      res$Q.resid <- res$df.Q.resid <- res$pval.Q.resid <- NA
      res$H.resid <- res$lower.H.resid <- res$upper.H.resid <- NA
      res$I2.resid <- res$lower.I2.resid <- res$upper.I2.resid <- NA
    }
    else {
      res$Q.w.random <- hcc$Q.resid
      res$df.Q.w.random <- hcc$df.Q.resid
      res$pval.Q.w.random <- hcc$pval.Q.resid
      ##
      res$tau2.resid <- hcc$tau2.resid
      res$lower.tau2.resid <- hcc$lower.tau2.resid
      res$upper.tau2.resid <- hcc$upper.tau2.resid
      ##
      res$tau.resid <- hcc$tau.resid
      res$lower.tau.resid <- hcc$lower.tau.resid
      res$upper.tau.resid <- hcc$upper.tau.resid
      res$sign.lower.tau.resid <- hcc$sign.lower.tau.resid
      res$sign.upper.tau.resid <- hcc$sign.upper.tau.resid
      ##
      res$Q.resid <- hcc$Q.resid
      res$df.Q.resid <- hcc$df.Q.resid
      res$pval.Q.resid <- hcc$pval.Q.resid
      ##
      res$H.resid <- hcc$H.resid
      res$lower.H.resid <- hcc$lower.H.resid
      res$upper.H.resid <- hcc$upper.H.resid
      ##
      res$I2.resid <- hcc$I2.resid
      res$lower.I2.resid <- hcc$lower.I2.resid
      res$upper.I2.resid <- hcc$upper.I2.resid
    }
    ##
    res$event.e.w <- NULL
    res$event.c.w <- NULL
    res$n.e.w <- NULL
    res$n.c.w <- NULL
    res$time.e.w <- NULL
    res$time.c.w <- NULL
  }
  ##
  ## Backward compatibility
  ##
  res$comb.fixed <- fixed
  res$comb.random <- random
  res$level.comb <- level.ma
  ##
  if (by) {
    res$byvar <- subgroup
    res$bylab <- subgroup.name
    res$print.byvar <- print.subgroup.name
    res$byseparator <- sep.subgroup
  }
  ##
  class(res) <- c(fun, "meta")

  
  res
}
