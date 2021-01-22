#' Meta-analysis of binary outcome data
#' 
#' @description
#' Calculation of fixed effect and random effects estimates (risk
#' ratio, odds ratio, risk difference, arcsine difference, or
#' diagnostic odds ratio) for meta-analyses with binary outcome
#' data. Mantel-Haenszel, inverse variance, Peto method, generalised
#' linear mixed model (GLMM), and sample size method are available for
#' pooling. For GLMMs, the \code{\link[metafor]{rma.glmm}} function
#' from R package \bold{metafor} (Viechtbauer, 2010) is called
#' internally.
#' 
#' @param event.e Number of events in experimental group or true
#'   positives in diagnostic study.
#' @param n.e Number of observations in experimental group or number
#'   of ill participants in diagnostic study.
#' @param event.c Number of events in control group or false positives
#'   in diagnostic study.
#' @param n.c Number of observations in control group or number of
#'   healthy participants in diagnostic study.
#' @param studlab An optional vector with study labels.
#' @param data An optional data frame containing the study
#'   information, i.e., event.e, n.e, event.c, and n.c.
#' @param subset An optional vector specifying a subset of studies to
#'   be used.
#' @param exclude An optional vector specifying studies to exclude
#'   from meta-analysis, however, to include in printouts and forest
#'   plots.
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"Inverse"},
#'   \code{"MH"}, \code{"Peto"}, \code{"GLMM"}, or \code{"SSW"}, can
#'   be abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"RR"}, \code{"OR"}, \code{"RD"}, \code{"ASD"}, or
#'   \code{"DOR"}) is to be used for pooling of studies, see Details.
#' @param incr Could be either a numerical value which is added to
#'   each cell frequency for studies with a zero cell count or the
#'   character string \code{"TACC"} which stands for treatment arm
#'   continuity correction, see Details.
#' @param allincr A logical indicating if \code{incr} is added to each
#'   cell frequency of all studies if at least one study has a zero
#'   cell count. If FALSE (default), \code{incr} is added only to each
#'   cell frequency of studies with a zero cell count.
#' @param addincr A logical indicating if \code{incr} is added to each
#'   cell frequency of all studies irrespective of zero cell counts.
#' @param allstudies A logical indicating if studies with zero or all
#'   events in both groups are to be included in the meta-analysis
#'   (applies only if \code{sm} is equal to \code{"RR"}, \code{"OR"},
#'   or \code{"DOR"}).
#' @param MH.exact A logical indicating if \code{incr} is not to be
#'   added to all cell frequencies for studies with a zero cell count
#'   to calculate the pooled estimate based on the Mantel-Haenszel
#'   method.
#' @param RR.Cochrane A logical indicating if 2*\code{incr} instead of
#'   1*\code{incr} is to be added to \code{n.e} and \code{n.c} in the
#'   calculation of the risk ratio (i.e., \code{sm="RR"}) for studies
#'   with a zero cell. This is used in RevMan 5, the program for
#'   preparing and maintaining Cochrane reviews.
#' @param Q.Cochrane A logical indicating if the Mantel-Haenszel
#'   estimate is used in the calculation of the heterogeneity
#'   statistic Q which is implemented in RevMan 5, the program for
#'   preparing and maintaining Cochrane reviews.
#' @param model.glmm A character string indicating which GLMM should
#'   be used.  One of \code{"UM.FS"}, \code{"UM.RS"}, \code{"CM.EL"},
#'   and \code{"CM.AL"}, see Details.
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
#'   funnel plot asymmetry is to be used. Either \code{"rank"},
#'   \code{"linreg"}, \code{"mm"}, \code{"count"}, \code{"score"}, or
#'   \code{"peters"}, can be abbreviated. See function
#'   \code{\link{metabias}.}
#' @param backtransf A logical indicating whether results for odds
#'   ratio (\code{sm="OR"}), risk ratio (\code{sm="RR"}), or
#'   diagnostic odds ratio (\code{sm="DOR"}) should be back
#'   transformed in printouts and plots. If TRUE (default), results
#'   will be presented as odds ratios and risk ratios; otherwise log
#'   odds ratios and log risk ratios will be shown.
#' @param pscale A numeric defining a scaling factor for printing of
#'   risk differences.
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
#' @param print.CMH A logical indicating whether result of the
#'   Cochran-Mantel-Haenszel test for overall effect should be
#'   printed.
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
#' with binary outcome data.
#' 
#' The following measures of treatment effect are available (Rücker et
#' al., 2009):
#' 
#' \itemize{
#' \item Risk ratio (\code{sm = "RR"})
#' \item Odds ratio (\code{sm = "OR"})
#' \item Risk difference (\code{sm = "RD"})
#' \item Arcsine difference (\code{sm = "ASD"})
#' \item Diagnostic Odds ratio (\code{sm = "DOR"})
#' }
#'
#' Note, mathematically, odds ratios and diagnostic odds ratios are
#' identical, however, the labels in printouts and figures differ.
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
#' Mantel-Haenszel method (Greenland & Robins, 1985; Robins et al.,
#' 1986) is used to calculate the fixed effect estimate; if
#' \code{method} is \code{"Inverse"}, inverse variance weighting is
#' used for pooling (Fleiss, 1993); if \code{method} is \code{"Peto"},
#' the Peto method is used for pooling (Yussuf et al., 1985); if
#' \code{method} is \code{"SSW"}, the sample size method is used for
#' pooling (Bakbergenuly et al., 2020).
#'
#' While the Mantel-Haenszel and Peto method are defined under the
#' fixed effect model, random effects variants based on these methods
#' are also implemented in \code{metabin}. Following RevMan 5, the
#' Mantel-Haenszel estimator is used in the calculation of the
#' between-study heterogeneity statistic Q which is used in the
#' DerSimonian-Laird estimator. Accordlingly, the results for the
#' random effects meta-analysis using the Mantel-Haenszel or inverse
#' variance method are typically very similar. For the Peto method,
#' Peto's log odds ratio, i.e. \code{(O-E) / V} and its standard error
#' \code{sqrt(1 / V)} with \code{O-E} and \code{V} denoting
#' "Observed minus Expected" and its variance, are utilised in the
#' random effects model. Accordingly, results of a random effects
#' model using \code{sm = "Peto"} can be different to results from a
#' random effects model using \code{sm = "MH"} or \code{sm =
#' "Inverse"}.
#' 
#' A distinctive and frequently overlooked advantage of binary
#' endpoints is that individual patient data (IPD) can be extracted
#' from a two-by-two table.  Accordingly, statistical methods for IPD,
#' i.e., logistic regression and generalised linear mixed models, can
#' be utilised in a meta-analysis of binary outcomes (Stijnen et al.,
#' 2010; Simmonds et al., 2016). These methods are available (argument
#' \code{method = "GLMM"}) for the odds ratio as summary measure by
#' calling the \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor} internally.
#'
#' Four different GLMMs are available for
#' meta-analysis with binary outcomes using argument \code{model.glmm}
#' (which corresponds to argument \code{model} in the
#' \code{\link[metafor]{rma.glmm}} function):
#' \tabular{cl}{
#' 1. \tab Logistic regression model with fixed study effects
#'  (default) \cr
#'  \tab (\code{model.glmm = "UM.FS"}, i.e., \bold{U}nconditional
#'  \bold{M}odel - \bold{F}ixed \bold{S}tudy effects) \cr
#' 2. \tab Mixed-effects logistic regression model with random study
#'  effects \cr
#'  \tab (\code{model.glmm = "UM.RS"}, i.e., \bold{U}nconditional
#'  \bold{M}odel - \bold{R}andom \bold{S}tudy effects) \cr
#' 3. \tab Generalised linear mixed model (conditional
#'  Hypergeometric-Normal) \cr
#'  \tab (\code{model.glmm = "CM.EL"}, i.e., \bold{C}onditional
#'  \bold{M}odel - \bold{E}xact \bold{L}ikelihood) \cr
#' 4. \tab Generalised linear mixed model (conditional
#'   Binomial-Normal) \cr
#'  \tab (\code{model.glmm = "CM.AL"}, i.e., \bold{C}onditional
#'   \bold{M}odel - \bold{A}pproximate \bold{L}ikelihood)
#' }
#'
#' Details on these four GLMMs as well as additional arguments which
#' can be provided using argument '\code{\dots}' in \code{metabin}
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
#' cell frequencies of these studies; if \code{incr} is \code{"TACC"}
#' a treatment arm continuity correction is used instead (Sweeting et
#' al., 2004; Diamond et al., 2007). For odds ratio and risk ratio,
#' treatment estimates and standard errors are only calculated for
#' studies with zero or all events in both groups if \code{allstudies}
#' is \code{TRUE}. This continuity correction is used both to
#' calculate individual study results with confidence limits and to
#' conduct meta-analysis based on the inverse variance method. For
#' Peto method and GLMMs no continuity correction is used. For the
#' Mantel-Haenszel method, by default (if \code{MH.exact} is FALSE),
#' \code{incr} is added to all cell frequencies of a study with a zero
#' cell count in the calculation of the pooled risk ratio or odds
#' ratio as well as the estimation of the variance of the pooled risk
#' difference, risk ratio or odds ratio. This approach is also used in
#' other software, e.g. RevMan 5 and the Stata procedure
#' metan. According to Fleiss (in Cooper & Hedges, 1994), there is no
#' need to add 0.5 to a cell frequency of zero to calculate the
#' Mantel-Haenszel estimate and he advocates the exact method
#' (\code{MH.exact} = TRUE). Note, estimates based on exact
#' Mantel-Haenszel method or GLMM are not defined if the number of
#' events is zero in all studies either in the experimental or control
#' group.
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
#' Wiksten et al., 2016). In such cases, an
#' \emph{ad hoc} variance correction has been proposed by utilising
#' the variance estimate from the classic random effects model with
#' the HK method (Knapp and Hartung, 2003; IQWiQ, 2020). An
#' alternative approach is to use the wider confidence interval of
#' classic fixed or random effects meta-analysis and the HK method
#' (Wiksten et al., 2016; Jackson et al., 2017).
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
#' An object of class \code{c("metabin", "meta")} with corresponding
#' \code{print}, \code{summary}, and \code{forest} functions. The
#' object is a list containing the following components:
#'
#' \item{event.e, n.e, event.c, n.c, studlab, exclude,}{As defined
#'   above.}
#' \item{sm, method, incr, allincr, addincr,}{As defined above.}
#' \item{allstudies, MH.exact, RR.Cochrane, Q.Cochrane, model.glmm,}{As
#'   defined above.}
#' \item{warn, level, level.comb, comb.fixed, comb.random,}{As defined
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
#' \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect,
#'   e.g., log risk ratio or risk difference, and standard error
#'   (fixed effect model).}
#' \item{lower.fixed, upper.fixed}{Lower and upper confidence interval
#'   limits (fixed effect model).}
#' \item{statistic.fixed, pval.fixed}{z-value and p-value for test of
#'   overall treatment effect (fixed effect model).}
#' \item{TE.random, seTE.random}{Estimated overall treatment effect,
#'   e.g., log risk ratio or risk difference, and standard error
#'   (random effects model).}
#' \item{lower.random, upper.random}{Lower and upper confidence
#'   interval limits (random effects model).}
#' \item{statistic.random, pval.random}{z-value or t-value and
#'   corresponding p-value for test of overall treatment effect
#'   (random effects model).}  \item{prediction, level.predict}{As
#'   defined above.}
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
#' \item{Q.CMH}{Cochran-Mantel-Haenszel test statistic for overall
#'   effect.}
#' \item{df.Q.CMH}{Degrees of freedom for Cochran-Mantel-Haenszel test
#'   statistic.}
#' \item{pval.Q.CMH}{P-value of Cochran-Mantel-Haenszel test.}
#' \item{incr.e, incr.c}{Increment added to cells in the experimental
#'   and control group, respectively.}
#' \item{sparse}{Logical flag indicating if any study included in
#'   meta-analysis has any zero cell frequencies.}
#' \item{doublezeros}{Logical flag indicating if any study has zero
#'   cell frequencies in both treatment groups.}
#' \item{df.hakn}{Degrees of freedom for test of treatment effect for
#'   Hartung-Knapp method (only if \code{hakn = TRUE}).}
#' \item{k.MH}{Number of studies combined in meta-analysis using
#'   Mantel-Haenszel method.}
#' \item{bylevs}{Levels of grouping variable - if \code{byvar} is not
#'   missing.}
#' \item{TE.fixed.w, seTE.fixed.w}{Estimated treatment effect and
#'   standard error in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}  \item{lower.fixed.w,
#'   upper.fixed.w}{Lower and upper confidence interval limits in
#'   subgroups (fixed effect model) - if \code{byvar} is not missing.}
#'
#' \item{statistic.fixed.w, pval.fixed.w}{z-value and p-value for test
#'   of treatment effect in subgroups (fixed effect model) - if
#'   \code{byvar} is not missing.}  \item{TE.random.w,
#'   seTE.random.w}{Estimated treatment effect and standard error in
#'   subgroups (random effects model) - if \code{byvar} is not
#'   missing.}
#' \item{lower.random.w, upper.random.w}{Lower and upper confidence
#'   interval limits in subgroups (random effects model) - if
#'   \code{byvar} is not missing.}
#' \item{statistic.random.w, pval.random.w}{z-value or t-value and
#'   corresponding p-value for test of treatment effect in subgroups
#'   (random effects model) - if \code{byvar} is not missing.}
#' \item{w.fixed.w, w.random.w}{Weight of subgroups (in fixed and
#'   random effects model) - if \code{byvar} is not missing.}
#'   \item{df.hakn.w}{Degrees of freedom for test of treatment effect
#'   for Hartung-Knapp method in subgroups - if \code{byvar} is not
#'   missing and \code{hakn = TRUE}.}  \item{event.e.w}{Number of
#'   events in experimental group in subgroups - if \code{byvar} is
#'   not missing.}  \item{n.e.w}{Number of observations in
#'   experimental group in subgroups - if \code{byvar} is not
#'   missing.}  \item{event.c.w}{Number of events in control group in
#'   subgroups - if \code{byvar} is not missing.}  \item{n.c.w}{Number
#'   of observations in control group in subgroups - if \code{byvar}
#'   is not missing.}  \item{k.w}{Number of studies combined within
#'   subgroups - if \code{byvar} is not missing.}
#'   \item{k.all.w}{Number of all studies in subgroups - if
#'   \code{byvar} is not missing.}  \item{Q.w.fixed}{Overall within
#'   subgroups heterogeneity statistic Q (based on fixed effect model)
#'   - if \code{byvar} is not missing.}  \item{Q.w.random}{Overall
#'   within subgroups heterogeneity statistic Q (based on random
#'   effects model) - if \code{byvar} is not missing (only calculated
#'   if argument \code{tau.common} is TRUE).}  \item{df.Q.w}{Degrees
#'   of freedom for test of overall within subgroups heterogeneity -
#'   if \code{byvar} is not missing.}  \item{pval.Q.w.fixed}{P-value
#'   of within subgroups heterogeneity statistic Q (based on fixed
#'   effect model) - if \code{byvar} is not missing.}
#'   \item{pval.Q.w.random}{P-value of within subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{byvar} is
#'   not missing.}  \item{Q.b.fixed}{Overall between subgroups
#'   heterogeneity statistic Q (based on fixed effect model) - if
#'   \code{byvar} is not missing.}  \item{Q.b.random}{Overall between
#'   subgroups heterogeneity statistic Q (based on random effects
#'   model) - if \code{byvar} is not missing.}  \item{df.Q.b}{Degrees
#'   of freedom for test of overall between subgroups heterogeneity -
#'   if \code{byvar} is not missing.}  \item{pval.Q.b.fixed}{P-value
#'   of between subgroups heterogeneity statistic Q (based on fixed
#'   effect model) - if \code{byvar} is not missing.}
#'   \item{pval.Q.b.random}{P-value of between subgroups heterogeneity
#'   statistic Q (based on random effects model) - if \code{byvar} is
#'   not missing.}  \item{tau.w}{Square-root of between-study variance
#'   within subgroups - if \code{byvar} is not missing.}
#'   \item{H.w}{Heterogeneity statistic H within subgroups - if
#'   \code{byvar} is not missing.}  \item{lower.H.w, upper.H.w}{Lower
#'   and upper confidence limit for heterogeneity statistic H within
#'   subgroups - if \code{byvar} is not missing.}
#'   \item{I2.w}{Heterogeneity statistic I\eqn{^2} within subgroups -
#'   if \code{byvar} is not missing.}  \item{lower.I2.w,
#'   upper.I2.w}{Lower and upper confidence limit for heterogeneity
#'   statistic I\eqn{^2} within subgroups - if \code{byvar} is not
#'   missing.}  \item{keepdata}{As defined above.}
#'   \item{data}{Original data (set) used in function call (if
#'   \code{keepdata = TRUE}).}  \item{subset}{Information on subset of
#'   original data used in meta-analysis (if \code{keepdata = TRUE}).}
#'   \item{.glmm.fixed}{GLMM object generated by call of
#'   \code{\link[metafor]{rma.glmm}} function (fixed effect model).}
#'   \item{.glmm.random}{GLMM object generated by call of
#'   \code{\link[metafor]{rma.glmm}} function (random effects model).}
#'   \item{call}{Function call.}  \item{version}{Version of R package
#'   \bold{meta} used to create object.}
#'   \item{version.metafor}{Version of R package \bold{metafor} used
#'   for GLMMs.}
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{update.meta}}, \code{\link{forest}},
#'   \code{\link{funnel}}, \code{\link{metabias}},
#'   \code{\link{metacont}}, \code{\link{metagen}},
#'   \code{\link{metareg}}, \code{\link{print.meta}}
#' 
#' @references
#' Bakbergenuly I, Hoaglin DC, Kulinskaya E (2020):
#' Methods for estimating between-study variance and overall
#' effect in meta-analysis of odds-ratios.
#' \emph{Research Synthesis Methods},
#' DOI: 10.1002/jrsm.1404
#' 
#' Cooper H & Hedges LV (1994):
#' \emph{The Handbook of Research Synthesis}.
#' Newbury Park, CA: Russell Sage Foundation
#' 
#' Diamond GA, Bax L, Kaul S (2007):
#' Uncertain Effects of Rosiglitazone on the Risk for Myocardial
#' Infarction and Cardiovascular Death.
#' \emph{Annals of Internal Medicine},
#' \bold{147}, 578--81
#' 
#' DerSimonian R & Laird N (1986):
#' Meta-analysis in clinical trials.
#' \emph{Controlled Clinical Trials},
#' \bold{7}, 177--88
#' 
#' Fleiss JL (1993):
#' The statistical basis of meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{2}, 121--45
#' 
#' Greenland S & Robins JM (1985):
#' Estimation of a common effect parameter from sparse follow-up data.
#' \emph{Biometrics},
#' \bold{41}, 55--68
#' 
#' Hartung J & Knapp G (2001):
#' A refined method for the meta-analysis of controlled clinical
#' trials with binary outcome.
#' \emph{Statistics in Medicine},
#' \bold{20}, 3875--89
#' 
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009):
#' A re-evaluation of random-effects meta-analysis.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \bold{172}, 137--59
#' 
#' IQWiG (2020):
#' General Methods: Version 6.0.
#' \url{https://www.iqwig.de/en/about-us/methods/methods-paper/}
#' 
#' Knapp G & Hartung J (2003):
#' Improved tests for a random effects meta-regression with a single
#' covariate.
#' \emph{Statistics in Medicine},
#' \bold{22}, 2693--710
#' 
#' \emph{Review Manager (RevMan)} [Computer program]. Version 5.4.
#' The Cochrane Collaboration, 2020
#' 
#' Paule RC & Mandel J (1982):
#' Consensus values and weighting factors.
#'\emph{Journal of Research of the National Bureau of Standards},
#' \bold{87}, 377--85
#' 
#' Pettigrew HM, Gart JJ, Thomas DG (1986):
#' The bias and higher cumulants of the logarithm of a binomial
#' variate.
#' \emph{Biometrika},
#' \bold{73}, 425--35
#'
#' Robins J, Breslow N, Greenland S (1986):
#' Estimators of the Mantel-Haenszel Variance Consistent in Both
#' Sparse Data and Large-Strata Limiting Models.
#' \emph{Biometrics},
#' \bold{42}, 311--23
#' 
#' Rücker G, Schwarzer G, Carpenter J, Olkin I (2009):
#' Why add anything to nothing? The arcsine difference as a measure of
#' treatment effect in meta-analysis with zero cells.
#' \emph{Statistics in Medicine},
#' \bold{28}, 721--38
#' 
#' Simmonds MC, Higgins JP (2016):
#' A general framework for the use of logistic regression models in
#' meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{25}, 2858--77
#' 
#' StataCorp. 2011.
#' \emph{Stata Statistical Software: Release 12}.
#' College Station, TX: StataCorp LP.
#' 
#' Stijnen T, Hamza TH, Ozdemir P (2010):
#' Random effects meta-analysis of event outcome in the framework of
#' the generalized linear mixed model with applications in sparse
#' data.
#' \emph{Statistics in Medicine},
#' \bold{29}, 3046--67
#' 
#' Sweeting MJ, Sutton AJ, Lambert PC (2004):
#' What to add to nothing? Use and avoidance of continuity corrections
#' in meta-analysis of sparse data.
#'\emph{Statistics in Medicine},
#' \bold{23}, 1351--75
#' 
#' Viechtbauer W (2010):
#' Conducting meta-analyses in R with the metafor package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#' 
#' Wiksten A, Rücker G, Schwarzer G (2016):
#' Hartung-Knapp method is not always conservative compared with
#' fixed-effect meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{35}, 2503--15
#'
#' Yusuf S, Peto R, Lewis J, Collins R, Sleight P (1985):
#' Beta blockade during and after myocardial infarction: An overview
#' of the randomized trials.
#' \emph{Progress in Cardiovascular Diseases},
#' \bold{27}, 335--71
#' 
#' @examples
#' # Calculate odds ratio and confidence interval for a single study
#' #
#' metabin(10, 20, 15, 20, sm = "OR")
#' 
#' # Different results (due to handling of studies with double zeros)
#' #
#' metabin(0, 10, 0, 10, sm = "OR")
#' metabin(0, 10, 0, 10, sm = "OR", allstudies = TRUE)
#' 
#' # Use subset of Olkin (1995) to conduct meta-analysis based on
#' # inverse variance method (with risk ratio as summary measure)
#' #
#' data(Olkin1995)
#' m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'               data = Olkin1995, subset = c(41, 47, 51, 59),
#'               method = "Inverse")
#' summary(m1)
#' 
#' # Use different subset of Olkin (1995)
#' #
#' m2 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'               data = Olkin1995, subset = year < 1970,
#'               method = "Inverse", studlab = author)
#' summary(m2)
#' forest(m2)
#' 
#' # Meta-analysis with odds ratio as summary measure
#' #
#' m3 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'               data = Olkin1995, subset = year < 1970,
#'               sm = "OR", method = "Inverse", studlab = author)
#' # Same meta-analysis result using 'update.meta' function
#' m3 <- update(m2, sm = "OR")
#' summary(m3)
#' 
#' # Meta-analysis based on Mantel-Haenszel method (with odds ratio as
#' # summary measure)
#' #
#' m4 <- update(m3, method = "MH")
#' summary(m4)
#' 
#' # Meta-analysis based on Peto method (only available for odds ratio
#' # as summary measure)
#' #
#' m5 <- update(m3, method = "Peto")
#' summary(m5)
#' 
#' \dontrun{
#' # Meta-analysis using generalised linear mixed models (only if R
#' # packages 'metafor' and 'lme4' are available)
#' #
#' if (suppressMessages(require(metafor, quietly = TRUE, warn = FALSE)) &
#'     require(lme4, quietly = TRUE)) {
#' 
#' # Logistic regression model with (k = 4) fixed study effects
#' # (default: model.glmm = "UM.FS")
#' #
#' m6 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'               data = Olkin1995, subset = year < 1970,
#'               method = "GLMM")
#' # Same results:
#' m6 <- update(m2, method = "GLMM")
#' summary(m6)
#' 
#' # Mixed-effects logistic regression model with random study effects
#' # (warning message printed due to argument 'nAGQ')
#' #
#' m7 <- update(m6, model.glmm = "UM.RS")
#' #
#' # Use additional argument 'nAGQ' for internal call of 'rma.glmm'
#' # function
#' #
#' m7 <- update(m6, model.glmm = "UM.RS", nAGQ = 1)
#' summary(m7)
#' 
#' # Generalised linear mixed model (conditional
#' # Hypergeometric-Normal) (R package 'BiasedUrn' must be available)
#' #
#' if (require(BiasedUrn, quietly = TRUE)) {
#'  m8 <- update(m6, model.glmm = "CM.EL")
#'  summary(m8)
#' }
#' 
#' # Generalised linear mixed model (conditional Binomial-Normal)
#' #
#' m9 <- update(m6, model.glmm = "CM.AL")
#' summary(m9)
#' 
#' # Logistic regression model with (k = 70) fixed study effects
#' # (about 18 seconds with Intel Core i7-3667U, 2.0GHz)
#' #
#' m10 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'                data = Olkin1995, method = "GLMM")
#' summary(m10)
#' 
#' # Mixed-effects logistic regression model with random study effects
#' # - about 50 seconds with Intel Core i7-3667U, 2.0GHz
#' # - several warning messages, e.g. "failure to converge, ..."
#' #
#' summary(update(m10, model.glmm = "UM.RS"))
#' 
#' # Conditional Hypergeometric-Normal GLMM
#' # - long computation time (about 12 minutes with Intel Core
#' #   i7-3667U, 2.0GHz)
#' # - estimation problems for this very large dataset:
#' #   * warning that Choleski factorization of Hessian failed
#' #   * confidence interval for treatment effect smaller in random
#' #     effects model compared to fixed effect model
#' #
#' if (require(BiasedUrn, quietly = TRUE)) {
#'  system.time(m11 <- update(m10, model.glmm = "CM.EL"))
#'  summary(m11)
#' }
#' 
#' # Generalised linear mixed model (conditional Binomial-Normal)
#' # (less than 1 second with Intel Core i7-3667U, 2.0GHz)
#' #
#' summary(update(m10, model.glmm = "CM.AL"))
#' }
#' }
#' 
#' @export metabin


metabin <- function(event.e, n.e, event.c, n.c, studlab,
                    ##
                    data = NULL, subset = NULL, exclude = NULL,
                    ##
                    method = ifelse(tau.common, "Inverse", gs("method")),
                    sm =
                      ifelse(!is.na(charmatch(tolower(method),
                                              c("peto", "glmm", "ssw"),
                                              nomatch = NA)),
                             "OR", gs("smbin")),
                    incr = gs("incr"), allincr = gs("allincr"),
                    addincr = gs("addincr"),
                    allstudies = gs("allstudies"),
                    MH.exact = gs("MH.exact"), RR.Cochrane = gs("RR.Cochrane"),
                    Q.Cochrane =
                      gs("Q.Cochrane") & method == "MH" & method.tau == "DL",
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
                    method.tau.ci = gs("method.tau.ci"),
                    tau.preset = NULL, TE.tau = NULL,
                    tau.common = gs("tau.common"),
                    ##
                    prediction = gs("prediction"),
                    level.predict = gs("level.predict"),
                    ##
                    method.bias = ifelse(sm == "OR", "score",
                                  ifelse(sm == "DOR", "deeks",
                                         gs("method.bias"))),
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
                    label.e = gs("label.e"), label.c = gs("label.c"),
                    label.left = gs("label.left"), label.right = gs("label.right"),
                    ##
                    byvar, bylab, print.byvar = gs("print.byvar"),
                    byseparator = gs("byseparator"),
                    ##
                    print.CMH = gs("print.CMH"),
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
  sm.metafor <- c("PHI", "YUQ", "YUY", "RTET",
                  "PBIT", "OR2D", "OR2DN", "OR2DL",
                  "MPRD", "MPRR", "MPOR", "MPORC", "MPPETO")
  ## sm <- setchar(sm, c(.settings$sm4bin, sm.metafor))
  sm <- setchar(sm, .settings$sm4bin)
  metafor <- sm %in% sm.metafor
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
  method.tau <- setchar(method.tau, c(.settings$meth4tau, "KD"))
  if (is.null(method.tau.ci))
    method.tau.ci <- if (method.tau == "DL") "J" else "QP"
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
  chknumeric(pscale, length = 1)
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
  ## Additional arguments / checks for metabin objects
  ##
  fun <- "metabin"
  ##
  chklogical(warn)
  if (sm != "RD" & pscale != 1) {
    if (warn)
      warning("Argument 'pscale' only considered for risk differences.")
    pscale <- 1
  }
  ##
  method <- setchar(method, .settings$meth4bin)
  if (metafor)
    method <- "Inverse"
  is.glmm <- method == "GLMM"
  ##
  chklogical(allincr)
  chklogical(addincr)
  chklogical(allstudies)
  chklogical(MH.exact)
  chklogical(RR.Cochrane)
  chklogical(Q.Cochrane)
  if (Q.Cochrane & (method != "MH" | method.tau != "DL")) {
    warning("Argument 'Q.Cochrane' only considered for ",
            "Mantel-Haenszel method in combination with ",
            "DerSimonian-Laird estimator.")
    Q.Cochrane <- FALSE
  }
  ##
  model.glmm <- setchar(model.glmm, c("UM.FS", "UM.RS", "CM.EL", "CM.AL"))
  if (is.glmm & model.glmm == "CM.EL")
    is.installed.package("BiasedUrn", fun, "model.glmm", " = \"CM.EL\"")
  ##
  chklogical(print.CMH)
  ##
  if (sm == "ASD") {
    method <- "Inverse"
    if (!missing(Q.Cochrane) && Q.Cochrane)
      warning("Argument 'Q.Cochrane' only considered for ",
              "Mantel-Haenszel method in combination with ",
              "DerSimonian-Laird estimator.")
    Q.Cochrane <- FALSE
  }
  ##
  if (sm != "OR") {
    if (method == "Peto")
      stop("Peto's method only possible with argument 'sm = \"OR\"'")
    else if (method == "SSW")
      stop("Sample size weighting only available with argument 'sm = \"OR\"'")
    else if (is.glmm)
      stop("Generalised linear mixed models only possible with ",
           "argument 'sm = \"OR\"'.")
  }
  ##
  if (is.glmm & method.tau != "ML")
    stop("Generalised linear mixed models only possible with ",
         "argument 'method.tau = \"ML\"'.")
  ##
  ## Check for deprecated arguments in '...'
  ##
  args  <- list(...)
  ## Check whether first argument is a list. In this case only use
  ## this list as input.
  if (length(args) > 0 && is.list(args[[1]]))
    args <- args[[1]]
  ##
  additional.arguments <- names(args)
  ##
  if (length(additional.arguments) > 0) {
    if ("RR.cochrane" %in% additional.arguments)
      if (!missing(RR.Cochrane))
        warning("Argument 'RR.cochrane' ignored as both arguments ",
                "'RR.Cochrane' and 'RR.cochrane' are provided.")
      else {
        RR.Cochrane <- args[["RR.cochrane"]]
        chklogical(RR.Cochrane)
      }
  }
  
  
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
  ## Catch 'event.e', 'n.e', 'event.c', and 'n.c' from data:
  ##
  event.e <- eval(mf[[match("event.e", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.e)
  k.All <- length(event.e)
  ##
  n.e <- eval(mf[[match("n.e", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.e)
  ##
  event.c <- eval(mf[[match("event.c", names(mf))]],
                  data, enclos = sys.frame(sys.parent()))
  chknull(event.c)
  ##
  n.c <- eval(mf[[match("n.c", names(mf))]],
              data, enclos = sys.frame(sys.parent()))
  chknull(n.c)
  ##
  ## Catch 'incr' from data:
  ##
  if (!missing(incr))
    incr <- eval(mf[[match("incr", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
  ##
  if (is.numeric(incr))
    chknumeric(incr, min = 0)
  else
    incr <- setchar(incr, "TACC",
                    "should be numeric or the character string \"TACC\"")
  ##
  if (metafor) {
    if (length(incr) > 1) {
      if (!missing(incr))
        warning("Increment of 0.5 used for effect measure '", sm, "'")
      incr <- 0.5
    }
    else if (incr == "TACC") {
      if (!missing(incr))
        warning("Increment of 0.5 used for effect measure '", sm, "'")
      incr <- 0.5
    }
  }
  ##
  ## Catch 'studlab', 'byvar', 'subset' and 'exclude' from data:
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
  chklength(n.e, k.All, fun)
  chklength(event.c, k.All, fun)
  chklength(n.c, k.All, fun)
  chklength(studlab, k.All, fun)
  ##
  if (length(incr) > 1)
    chklength(incr, k.All, fun)
  ##
  if (by)
    chklength(byvar, k.All, fun)
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
    if (warn)
      warning("Value for argument 'tau.common' set to FALSE as argument 'byvar' is missing.")
    tau.common <- FALSE
  }
  if (by & !tau.common & !is.null(tau.preset)) {
    if (warn)
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
      data <- data.frame(.event.e = event.e)
    else
      data$.event.e <- event.e
    ##
    data$.n.e <- n.e
    data$.event.c <- event.c
    data$.n.c <- n.c
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
  }


  ##
  ##
  ## (6) Use subset for analysis
  ##
  ##
  if (!missing.subset) {
    event.e <- event.e[subset]
    n.e <- n.e[subset]
    event.c <- event.c[subset]
    n.c <- n.c[subset]
    studlab <- studlab[subset]
    ##
    exclude <- exclude[subset]
    ##
    if (length(incr) > 1)
      incr <- incr[subset]
    ##
    if (by)
      byvar <- byvar[subset]
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
  chknumeric(event.e)
  chknumeric(n.e)
  chknumeric(event.c)
  chknumeric(n.c)
  ##
  ## Recode integer as numeric:
  ##
  event.e <- int2num(event.e)
  n.e     <- int2num(n.e)
  event.c <- int2num(event.c)
  n.c     <- int2num(n.c)
  ##
  if (by) {
    chkmiss(byvar)
    byvar.name <- byvarname(mf[[match("byvar", names(mf))]])
    bylab <- if (!missing(bylab) && !is.null(bylab)) bylab else byvar.name
  }
  
  
  ##
  ##
  ## (7) Calculate results for individual studies
  ##
  ##
  ## Include non-informative studies?
  ## (i.e. studies with either zero or all events in both groups)
  ##
  if (sm == "RD" | sm == "ASD" | metafor)
    incl <- rep(1, k.all)
  else {
    allevents <- event.c == n.c & event.e == n.e
    if (allstudies)
      incl <- rep(1, k.all)
    else {
      if (sm %in% c("OR", "DOR"))
        incl <- ifelse((event.c == 0   & event.e == 0) |
                       (event.c == n.c & event.e == n.e), NA, 1)
      if (sm == "RR")
        incl <- ifelse((event.c == 0 & event.e == 0), NA, 1)
    }
  }
  ##
  ## Exclude studies from meta-analysis:
  ##
  sel1 <- event.e > n.e
  sel2 <- event.c > n.c
  if ((any(sel1, na.rm = TRUE)) & warn)
    warning("Studies with event.e > n.e get no weight in meta-analysis.")
  if ((any(sel2, na.rm = TRUE)) & warn)
    warning("Studies with event.c > n.c get no weight in meta-analysis.")
  incl[sel1 | sel2] <- NA
  ##
  sel3 <- n.e <= 0 | n.c <= 0
  if ((any(sel3, na.rm = TRUE)) & warn)
    warning("Studies with non-positive values for n.e and / or n.c ",
            "get no weight in meta-analysis.")
  incl[sel3] <- NA
  ##
  sel4 <- event.e < 0 | event.c < 0
  if ((any(sel4, na.rm = TRUE)) & warn)
    warning("Studies with negative values for event.e and / or event.c ",
            "get no weight in meta-analysis.")
  incl[sel4] <- NA
  ##
  ## Sparse computation
  ##
  sel <- switch(sm,
                OR = ((n.e - event.e) == 0 | event.e == 0 |
                      (n.c - event.c) == 0 | event.c == 0),
                RD = ((n.e - event.e) == 0 | event.e == 0 |
                      (n.c - event.c) == 0 | event.c == 0),
                RR = ((n.e - event.e) == 0 | event.e == 0 |
                      (n.c - event.c) == 0 | event.c == 0),
                ASD = rep(FALSE, length(event.e)),
                DOR = ((n.e - event.e) == 0 | event.e == 0 |
                       (n.c - event.c) == 0 | event.c == 0))
  ##
  sel[is.na(incl)] <- FALSE
  ##
  sparse <- any(sel, na.rm = TRUE)
  ##
  ## Check for studies with zero cell frequencies in both groups
  ##
  doublezeros <- FALSE
  if (sparse & sm %in% c("RR", "OR") & !(method %in% c("Peto", "GLMM"))) {
    sel.doublezeros <- switch(sm,
                              OR = (event.e == 0   & event.c ==   0) |
                                   (event.c == n.c & event.e == n.e),
                              RR = (event.c == 0 & event.e == 0))
    if (any(sel.doublezeros, na.rm = TRUE))
      doublezeros <- TRUE
  }
  ##
  ## No need to add anything to cell counts for
  ##  (i)  arcsine difference as summary measure
  ##  (ii) Peto method or GLMM
  ##
  if (sm == "ASD" | method %in% c("Peto", "GLMM")) {
    if ((!missing(incr) & any(incr != 0)) |
        (!missing(allincr) & allincr ) |
        (!missing(addincr) & addincr) |
        (!missing(allstudies) & allstudies)
        )
      if (sm == "ASD") {
        if ((sparse | addincr) & warn) {
          warning("Note, no continuity correction considered ",
                  "for arcsine difference (sm = \"ASD\").")
        }
      }
      else if (method == "Peto") {
        if ((sparse | addincr) & warn)
          warning("Note, no continuity correction considered ",
                  "for method = \"Peto\".")
      }
      else if (is.glmm) {
        if ((sparse | addincr) & warn)
          warning("Note, for method = \"GLMM\", continuity correction ",
                  "only used to calculate individual study results.")
      }
  }
  ##
  ## Define continuity correction
  ##
  if (addincr) {
    ##
    if (is.numeric(incr)) {
      incr.e <- if (length(incr) == 1) rep(incr, k.all) else incr
      incr.c <- if (length(incr) == 1) rep(incr, k.all) else incr
    }
    else {
      if (all(incr == "TACC")) {
        ##
        ## Treatment arm continuity correction:
        ##
        incr.e <- n.e / (n.e + n.c)
        incr.c <- n.c / (n.e + n.c)
      }
    }
  }
  else {
    if (sparse) {
      if (allincr) {
        ##
        if (is.numeric(incr)) {
          incr.e <- if (length(incr) == 1) rep(incr, k.all) else incr
          incr.c <- if (length(incr) == 1) rep(incr, k.all) else incr
        }
        else {
          if (all(incr == "TACC")) {
            ##
            ## Treatment arm continuity correction:
            ##
            incr.e <- n.e / (n.e + n.c)
            incr.c <- n.c / (n.e + n.c)
          }
        }
      }
      else {
        ##
        ## Bradburn, Deeks, Altman, Stata-procedure "metan":
        ## & SAS PROC FREQ (for method = "Inverse")
        ##
        if (is.numeric(incr)) {
          incr.e <- incr * sel
          incr.c <- incr * sel
        }
        else {
          if (all(incr == "TACC")) {
            ##
            ## Treatment arm continuity correction:
            ##
            incr.e <- n.e / (n.e + n.c) * sel
            incr.c <- n.c / (n.e + n.c) * sel
          }
        }
      }
    }
    else {
      incr.e <- rep(0, k.all)
      incr.c <- rep(0, k.all)
    }
  }
  ##
  ## No continuity correction for Peto method
  ##
  if (method == "Peto") {
    incr <- 0
    incr.e <- rep(0, k.all)
    incr.c <- rep(0, k.all)
  }
  ##
  n11 <- event.e * incl
  n21 <- event.c * incl
  n1. <- n.e * incl
  n2. <- n.c * incl
  ##
  n.. <- n1. + n2.
  n12 <- n1. - n11
  n22 <- n2. - n21
  n.1 <- n11 + n21
  n.2 <- n12 + n22
  ##
  Q.CMH <- (sum((n11 - n1. * n.1 / n..)[!exclude], na.rm = TRUE)^2 /
              sum((n1. * n2. * n.1 * n.2 / n..^3)[!exclude], na.rm = TRUE))
  ##
  p.e <- (n11 + incr.e) / (n1. + 2 * incr.e)
  p.c <- (n21 + incr.c) / (n2. + 2 * incr.c)
  ##
  ## Estimation of treatment effects in individual studies
  ##
  if (sm %in% c("OR", "DOR")) {
    if (method %in% c("MH", "Inverse", "GLMM", "SSW")) {
      ##
      ## Cooper & Hedges (1994), p. 251-2
      ##
      TE <- log(((n11 + incr.e) * (n22 + incr.c)) /
                ((n12 + incr.e) * (n21 + incr.c)))
      seTE <- sqrt((1 / (n11 + incr.e) + 1 / (n12 + incr.e) +
                    1 / (n21 + incr.c) + 1 / (n22 + incr.c)))
    }
    else if (method == "Peto") {
      ##
      ## Cooper & Hedges (1994), p. 252
      ##
      O <- n11
      E <- n1. * n.1 / n..
      V <- n1. * n2. * n.1 * n.2 / ((n.. - 1) * n..^2)
      ##
      TE <- (O - E) / V
      seTE <- sqrt(1 / V)
    }
  }
  else if (sm == "RR") {
    ##
    ## Cooper & Hedges (1994), p. 247-8
    ##
    if (!RR.Cochrane) {
      TE <- log(((n11 + incr.e) / (n1. + incr.e)) /
                  ((n21 + incr.c) / (n2. + incr.c)))
      ##
      ## Hartung & Knapp (2001), Stat Med, equation (18)
      ##
      seTE <- sqrt((1 / (n11 + incr.e * (!allevents)) - 1 / (n1. + incr.e) +
                      1 / (n21 + incr.c * (!allevents)) - 1 / (n2. + incr.c)))
    }
    else {
      TE <- log(((n11 + incr.e) / (n1. + 2 * incr.e)) /
                  ((n21 + incr.c) / (n2. + 2 * incr.c)))
      seTE <- sqrt((1 / (n11 + incr.e) - 1 / (n1. + 2 * incr.e) +
                      1 / (n21 + incr.c) - 1 / (n2. + 2 * incr.c)))
    }
  }
  else if (sm == "RD") {
    ##
    ## Cooper & Hedges (1994), p. 246-7
    ##
    TE <- n11 / n1. - n21 / n2.
    seTE <- sqrt((n11 + incr.e) * (n12 + incr.e) / (n1. + 2 * incr.e)^3 +
                   (n21 + incr.c) * (n22 + incr.c) / (n2. + 2 * incr.c)^3)
  }
  else if (sm == "ASD") {
    ##
    ## Ruecker et al. (2009)
    ##
    TE <- asin(sqrt(n11 / n1.)) - asin(sqrt(n21 / n2.))
    seTE <- sqrt(0.25 * (1 / n1. + 1 / n2.))
  }
  else if (metafor) {
    ##
    ## Other effect measures calculated in R package metafor
    ##
    if (addincr)
      to <- "all"
    else if (allincr)
      to <- "if0all"
    else
      to <- "only0"
    ##
    tmp <- escalc(measure = sm,
                  ai = n11, bi = n12, ci = n21, di = n22,
                  add = incr, to = to, drop00 = !allstudies)
    TE <- tmp$yi
    seTE <- sqrt(tmp$vi)
  }
  
  
  ##
  ##
  ## (8) Do meta-analysis
  ##
  ##
  k <- sum(!is.na(event.e[!exclude]) & !is.na(event.c[!exclude]) &
           !is.na(n.e[!exclude]) & !is.na(n.c[!exclude]))
  ##
  if (k == 1 & hakn)
    hakn <- FALSE
  ##
  if (all(incr.e == 0) & all(incr.c == 0) & method == "MH" & MH.exact == FALSE)
    MH.exact <- TRUE
  ##
  if (method == "MH") {
    incr.e <- incr.e * (!MH.exact)
    incr.c <- incr.c * (!MH.exact)
    ##
    if (sm %in% c("OR", "DOR")) {
      ##
      ## Cooper & Hedges (1994), p. 253-5 (MH.exact == TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## und RevMan 3.1 (MH.exact == FALSE)
      ##
      A <- (n11 + incr.e) * (n22 + incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      B <- (n11 + incr.e + n22 + incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      C <- (n12 + incr.e) * (n21 + incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      D <- (n12 + incr.e + n21 + incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      ##
      A[exclude] <- B[exclude] <- C[exclude] <- D[exclude] <- 0
      ##
      ## Cooper & Hedges (1994), p. 265-6
      ##
      w.fixed <- C
      TE.fixed <- log(sum(A, na.rm = TRUE) / sum(C, na.rm = TRUE))
      seTE.fixed <- sqrt((1 / (2 * sum(A, na.rm = TRUE)^2)  *
                            (sum(A * B, na.rm = TRUE) +
                               exp(TE.fixed) * (sum(B * C, na.rm = TRUE) +
                                                  sum(A * D, na.rm = TRUE)) +
                                                    exp(TE.fixed)^2 * sum(C * D, na.rm = TRUE))))
    }
    else if (sm == "RR") {
      ##
      ## Greenland, Robins (1985) (MH.exact == TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## (MH.exact == FALSE)
      ##
      D <- ((n1. + 2 * incr.e) * (n2. + 2 * incr.c) * (n.1 + incr.e + incr.c) -
              (n11 + incr.e) * (n21 + incr.c) * (n.. + 2 * incr.e + 2 * incr.c)) /
                (n.. + 2 * incr.e + 2 * incr.c)^2
      R <- (n11 + incr.e) * (n2. + 2 * incr.c) / (n.. + 2 * incr.e + 2 * incr.c)
      S <- (n21 + incr.c) * (n1. + 2 * incr.e) / (n.. + 2 * incr.e + 2 * incr.c)
      ##
      D[exclude] <- R[exclude] <- S[exclude] <- 0
      ##
      w.fixed <- S
      TE.fixed <- log(sum(R, na.rm = TRUE) / sum(S, na.rm = TRUE))
      seTE.fixed <- sqrt(sum(D, na.rm = TRUE) / (sum(R, na.rm = TRUE) *
                                                   sum(S, na.rm = TRUE)))
    }
    else if (sm == "RD") {
      ##
      ## Jon Deeks (1999) (MH.exact == TRUE)
      ##
      ## Bradburn, Deeks, Altman, Stata-procedure "metan"
      ## und RevMan 3.1 (MH.exact == FALSE)
      ##
      R <- ((n11 + incr.e) * (n12 + incr.e) * (n2. + 2 * incr.c)^3 +
              (n21 + incr.c) * (n22 + incr.c) * (n1. + 2 * incr.e)^3) /
                ((n1. + 2 * incr.e) * (n2. + 2 * incr.c) * (n.. + 2 * incr.e + 2 * incr.c)^2)
      S <- (n1. + 2 * incr.e) * (n2. + 2 * incr.c) /
        (n.. + 2 * incr.e + 2 * incr.c)
      ##
      R[exclude] <- S[exclude] <- 0
      ##
      w.fixed <- S
      TE.fixed <- weighted.mean(TE, w.fixed, na.rm = TRUE)
      seTE.fixed <- sqrt(sum(R, na.rm = TRUE) / sum(S, na.rm = TRUE)^2)
    }
    ##
    w.fixed[is.na(w.fixed)] <- 0
  }
  else if (method == "Peto") {
    w.fixed <- 1 / seTE^2
    w.fixed[exclude] <- 0
    TE.fixed   <- weighted.mean(TE, w.fixed, na.rm = TRUE)
    seTE.fixed <- sqrt(1 / sum(w.fixed, na.rm = TRUE))
    ##
    w.fixed[is.na(w.fixed)] <- 0
  }
  else if (is.glmm) {
    zero.all <-
      (sum(event.e[!exclude], na.rm = TRUE) == 0 &
       sum(event.c[!exclude], na.rm = TRUE) == 0) |
      (!any(event.e[!exclude] != n.e[!exclude]) |
       !any(event.c[!exclude] != n.c[!exclude]))
    ##
    if (!zero.all)
      glmm.fixed <- rma.glmm(ai = event.e[!exclude], n1i = n.e[!exclude],
                             ci = event.c[!exclude], n2i = n.c[!exclude],
                             method = "FE", test = ifelse(hakn, "t", "z"),
                             level = 100 * level.comb,
                             measure = "OR", model = model.glmm,
                             control = control,
                             ...)
    else
      glmm.fixed <- list(b = NA, se = NA,
                         QE.Wld = NA, QE.df = NA, QE.LRT = NA,
                         tau2 = NA, se.tau2 = NA)
    ##
    TE.fixed   <- as.numeric(glmm.fixed$b)
    seTE.fixed <- as.numeric(glmm.fixed$se)
    ##
    w.fixed <- rep(NA, length(event.e))
  }
  else if (method == "SSW") {
    w.fixed <- n.e * n.c / (n.e + n.c)
    w.fixed[exclude] <- 0
    TE.fixed <- weighted.mean(TE, w.fixed, na.rm = TRUE)
    seTE.fixed <- sqrt(sum(w.fixed^2 * seTE^2, na.rm = TRUE) /
                       sum(w.fixed, na.rm = TRUE)^2)
    ##
    w.fixed[is.na(w.fixed)] <- 0
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
               TE.tau = if (Q.Cochrane) TE.fixed else TE.tau,
               tau.common = FALSE,
               ##
               prediction = prediction,
               level.predict = level.predict,
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
               label.e = label.e, label.c = label.c,
               label.left = label.left, label.right = label.right,
               ##
               keepdata = FALSE,
               warn = warn,
               ##
               control = control)
  ##
  if (method == "SSW") {
    w.random <- n.e * n.c / (n.e + n.c)
    w.random[exclude] <- 0
    TE.random <- weighted.mean(TE, w.random, na.rm = TRUE)
    seTE.random <- sqrt(sum(w.random^2 * (seTE^2 + m$tau^2), na.rm = TRUE) /
                        sum(w.random, na.rm = TRUE)^2)
    ##
    w.random[is.na(w.random)] <- 0
  }
  ##
  if (by & tau.common & !is.glmm) {
    ## Estimate common tau-squared across subgroups
    hcc <- hetcalc(TE, seTE, method.tau, "",
                   if (Q.Cochrane & method == "MH") TE.fixed else TE.tau,
                   level.comb, byvar, control)
  }
  
  
  ##
  ##
  ## (9) Generate R object
  ##
  ##
  res <- list(event.e = event.e, n.e = n.e,
              event.c = event.c, n.c = n.c,
              method = method,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              sparse = sparse,
              allincr = allincr, addincr = addincr,
              allstudies = allstudies,
              doublezeros = doublezeros,
              MH.exact = MH.exact, RR.Cochrane = RR.Cochrane,
              Q.Cochrane = Q.Cochrane,
              Q.CMH = Q.CMH, df.Q.CMH = 1, pval.Q.CMH = pvalQ(Q.CMH, 1),
              print.CMH = print.CMH,
              incr.e = incr.e, incr.c = incr.c,
              k.MH = if (method == "MH") sum(w.fixed > 0) else NA,
              k.all = k.all)
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
  res$TE.tau <- TE.tau
  ##
  res$pscale <- pscale
  ##
  res$call <- match.call()
  ##
  if (method %in% c("MH", "Peto", "GLMM", "SSW")) {
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
    res$zval.fixed <- ci.f$statistic
  }
  ##
  if (is.glmm) {
    ##
    if (sum(!exclude) > 1 & !zero.all)
      glmm.random <- rma.glmm(ai = event.e[!exclude], n1i = n.e[!exclude],
                              ci = event.c[!exclude], n2i = n.c[!exclude],
                              method = method.tau,
                              test = ifelse(hakn, "t", "z"),
                              level = 100 * level.comb,
                                measure = "OR", model = model.glmm,
                              control = control,
                              ...)
    else {
      ##
      ## Fallback to fixed effect model due to small number of studies
      ## or zero or all events in all studies
      ##
      glmm.random <- glmm.fixed
    }
    ##
    TE.random   <- as.numeric(glmm.random$b)
    seTE.random <- as.numeric(glmm.random$se)
    ##
    ci.r <- ci(TE.random, seTE.random, level = level.comb,
               df = if (hakn) k - 1)
    ##
    res$w.random <- rep(NA, length(event.e))
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
    res$model.glmm <- model.glmm
    ##
    res$Q      <- glmm.random$QE.Wld
    res$df.Q   <- glmm.random$QE.df
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
    H <- calcH(res$Q, res$df.Q, level.comb)
    res$H <- H$TE
    res$lower.H <- H$lower
    res$upper.H <- H$upper
    ##
    I2 <- isquared(res$Q, res$df.Q, level.comb)
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
      n.by <- length(unique(byvar[!exclude]))
      if (n.by > 1)
        byvar.glmm <- factor(byvar[!exclude], bylevs(byvar[!exclude]))
      ##
      glmm.random.by <-
        try(suppressWarnings(rma.glmm(ai = event.e[!exclude],
                                      n1i = n.e[!exclude],
                                      ci = event.c[!exclude],
                                      n2i = n.c[!exclude],
                                      mods =
                                        if (n.by > 1) ~ byvar.glmm else NULL,
                                      method = method.tau,
                                      test = ifelse(hakn, "t", "z"),
                                      level = 100 * level.comb,
                                      measure = "OR", model = model.glmm,
                                      control = control,
                                      ...)),
            silent = TRUE)
      ##
      if ("try-error" %in% class(glmm.random.by))
        if (grepl(paste0("Number of parameters to be estimated is ",
                         "larger than the number of observations"),
                  glmm.random.by)) {
          glmm.random.by <-
            suppressWarnings(rma.glmm(ai = event.e[!exclude],
                                      n1i = n.e[!exclude],
                                      ci = event.c[!exclude],
                                      n2i = n.c[!exclude],
                                      mods =
                                        if (n.by > 1) ~ byvar.glmm else NULL,
                                      method = "FE",
                                      test = ifelse(hakn, "t", "z"),
                                      level = 100 * level.comb,
                                      measure = "OR", model = model.glmm,
                                      control = control,
                                      ...))
        }
        else
          stop(glmm.random.by)
      ##
      Q.r <- glmm.random.by$QE.Wld
      df.Q.r <- glmm.random.by$k - glmm.random.by$p
      ##
      H.r  <- calcH(Q.r, df.Q.r, level.comb)
      I2.r <- isquared(Q.r, df.Q.r, level.comb)
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
                  Q.resid = NA,
                  df.Q.resid = NA,
                  pval.Q.resid = NA,
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
  else if (method == "SSW") {
    if (hakn)
      ci.r <- ci(TE.random, seTE.random, level = level.comb, df = m$k - 1)
    else
      ci.r <- ci(TE.random, seTE.random, level = level.comb)
    ##
    res$TE.random <- TE.random
    res$seTE.random <- seTE.random
    res$w.random <- w.random
    res$lower.random <- ci.r$lower
    res$upper.random <- ci.r$upper
    res$statistic.random <- ci.r$statistic
    res$pval.random <- ci.r$p      
    res$zval.random <- ci.r$statistic
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
    if (!tau.common)
      res <- c(res, subgroup(res))
    else if (!is.null(tau.preset))
      res <- c(res, subgroup(res, tau.preset))
    else {
      if (is.glmm)
        res <- c(res, subgroup(res, NULL,
                               factor(res$byvar, bylevs(res$byvar)), ...))
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
      res$Q.w.random <- hcc$Q.resid
      res$df.Q.w.random <- hcc$df.Q.resid
      res$pval.Q.w.random <- hcc$pval.Q.resid
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
    res$event.w <- NULL
    res$n.w <- NULL
    res$time.e.w <- NULL
    res$time.c.w <- NULL
  }
  ##
  class(res) <- c(fun, "meta")


  res
}
