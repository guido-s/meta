#' Meta-analysis of binary outcome data
#' 
#' @description
#' Calculation of common effect and random effects estimates (risk
#' ratio, odds ratio, risk difference, arcsine difference, or
#' diagnostic odds ratio) for meta-analyses with binary outcome
#' data. Mantel-Haenszel, inverse variance, Peto method, generalised
#' linear mixed model (GLMM), logistic regression with penalised likelihood
#' and sample size method are available for pooling. For GLMMs,
#' the \code{\link[metafor]{rma.glmm}} function from R package \bold{metafor}
#' (Viechtbauer, 2010) is called internally. For penalised logistic regression,
#' R package \bold{brglm2} must be available.
#' 
#' @param event.e Number of events in experimental group, or true
#'   positives in diagnostic study, or an R object
#'   created with \code{\link{pairwise}}.
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
#' @param cluster An optional vector specifying which estimates come
#'   from the same cluster resulting in the use of a three-level
#'   meta-analysis model.
#' @param rho Assumed correlation of estimates within a cluster.
#' @param weights A single numeric or vector with user-specified weights.
#' @param weights.common User-specified weights (common effect model).
#' @param weights.random User-specified weights (random effects model).
#' @param method A character string indicating which method is to be
#'   used for pooling of studies. One of \code{"Inverse"},
#'   \code{"MH"}, \code{"Peto"}, \code{"GLMM"}, \code{"LRP"}, or \code{"SSW"},
#'   can be abbreviated.
#' @param sm A character string indicating which summary measure
#'   (\code{"RR"}, \code{"OR"}, \code{"RD"}, \code{"ASD"},
#'   \code{"DOR"}, or \code{"VE"}) is to be used for pooling of
#'   studies, see Details.
#' @param incr Could be either a numerical value which is added to
#'   cell frequencies for studies with a zero cell count, the
#'   character string \code{"TACC"} which stands for treatment arm
#'   continuity correction, or a numeric vector with the continuity
#'   correction for each study, see Details.
#' @param method.incr A character string indicating which continuity
#'   correction method should be used (\code{"only0"},
#'   \code{"if0all"}, \code{"all"}, or \code{"user"}), see Details.
#' @param allstudies A logical indicating if studies with zero or all
#'   events in both groups are to be included in the meta-analysis
#'   (applies only if \code{sm} is equal to \code{"RR"}, \code{"OR"},
#'   or \code{"DOR"}).
#' @param incr.e Continuity correction in experimental group, see Details.
#' @param incr.c Continuity correction in control group, see Details.
#' @param MH.exact A logical indicating if \code{incr} is not to be
#'   added to cell frequencies for studies with a zero cell count to
#'   calculate the pooled estimate based on the Mantel-Haenszel
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
#' @param common A logical indicating whether a common effect
#'   meta-analysis should be conducted.
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
#' @param method.tau A character string indicating which method is
#'   used to estimate the between-study variance \eqn{\tau^2} and its
#'   square root \eqn{\tau} (see \code{\link{meta-package}}).
#' @param method.tau.ci A character string indicating which method is
#'   used to estimate the confidence interval of \eqn{\tau^2} and
#'   \eqn{\tau} (see \code{\link{meta-package}}).
#' @param level.hetstat The level used to calculate confidence intervals
#'   for heterogeneity statistics.
#' @param tau.preset Prespecified value for the square root of the
#'   between-study variance \eqn{\tau^2}.
#' @param TE.tau Overall treatment effect used to estimate the
#'   between-study variance tau-squared.
#' @param tau.common A logical indicating whether tau-squared should
#'   be the same across subgroups.
#' @param detail.tau Detail on between-study variance estimate.
#' @param method.I2 A character string indicating which method is
#'   used to estimate the heterogeneity statistic I\eqn{^2}. Either
#'   \code{"Q"} or \code{"tau2"}, can be abbreviated
#'   (see \code{\link{meta-package}}).
#' @param level.ma The level used to calculate confidence intervals
#'   for meta-analysis estimates.
#' @param method.common.ci A character string indicating which method
#'   is used to calculate confidence interval and test statistic for
#'   common effect estimate (see \code{\link{meta-package}}).
#' @param method.random.ci A character string indicating which method
#'   is used to calculate confidence interval and test statistic for
#'   random effects estimate (see \code{\link{meta-package}}).
#' @param adhoc.hakn.ci A character string indicating whether an
#'   \emph{ad hoc} variance correction should be applied in the case
#'   of an arbitrarily small Hartung-Knapp variance estimate (see
#'   \code{\link{meta-package}}).
#' @param level.predict The level used to calculate prediction
#'   interval for a new study.
#' @param method.predict A character string indicating which method is
#'   used to calculate a prediction interval (see
#'   \code{\link{meta-package}}).
#' @param adhoc.hakn.pi A character string indicating whether an
#'   \emph{ad hoc} variance correction should be applied for
#'   prediction interval (see \code{\link{meta-package}}).
#' @param seed.predict A numeric value used as seed to calculate
#'   bootstrap prediction interval (see \code{\link{meta-package}}).
#' @param method.bias A character string indicating which test for
#'   funnel plot asymmetry is to be used. Either \code{"Begg"},
#'   \code{"Egger"}, \code{"Thompson"}, \code{"Schwarzer"},
#'   \code{"Harbord"}, \code{"Peters"}, or \code{"Deeks"}, can be
#'   abbreviated. See function \code{\link{metabias}.}
#' @param backtransf A logical indicating whether results for odds
#'   ratio (\code{sm="OR"}), risk ratio (\code{sm="RR"}), or
#'   diagnostic odds ratio (\code{sm="DOR"}) should be back
#'   transformed in printouts and plots. If TRUE (default), results
#'   will be presented as odds ratios and risk ratios; otherwise log
#'   odds ratios and log risk ratios will be shown.
#' @param pscale A numeric defining a scaling factor for printing of
#'   risk differences.
#' @param text.common A character string used in printouts and forest
#'   plot to label the pooled common effect estimate.
#' @param text.random A character string used in printouts and forest
#'   plot to label the pooled random effects estimate.
#' @param text.predict A character string used in printouts and forest
#'   plot to label the prediction interval.
#' @param text.w.common A character string used to label weights of
#'   common effect model.
#' @param text.w.random A character string used to label weights of
#'   random effects model.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' @param label.e Label for experimental group.
#' @param label.c Label for control group.
#' @param label.left Graph label on left side of null effect in forest plot.
#' @param label.right Graph label on right side of null effect in forest plot.
#' @param col.label.left The colour of the graph label on the left side of
#'   the null effect.
#' @param col.label.right The colour of the graph label on the right side of
#'   the null effect.
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
#' @param prediction.subgroup A logical indicating whether prediction
#'   intervals should be printed for subgroups.
#' @param seed.predict.subgroup A numeric vector providing seeds to
#'   calculate bootstrap prediction intervals within subgroups. Must
#'   be of same length as the number of subgroups.
#' @param byvar Deprecated argument (replaced by 'subgroup').
#' @param hakn Deprecated argument (replaced by 'method.random.ci').
#' @param adhoc.hakn Deprecated argument (replaced by
#'   'adhoc.hakn.ci').
#' @param print.CMH A logical indicating whether result of the
#'   Cochran-Mantel-Haenszel test for overall effect should be
#'   printed.
#' @param keepdata A logical indicating whether original data (set)
#'   should be kept in meta object.
#' @param warn A logical indicating whether warnings should be printed
#'   (e.g., if \code{incr} is added to studies with zero cell
#'   frequencies or if estimation problems exist in fitting a GLMM).
#' @param warn.deprecated A logical indicating whether warnings should
#'   be printed if deprecated arguments are used.
#' @param control An optional list to control the iterative process to
#'   estimate the between-study variance \eqn{\tau^2}. This argument
#'   is passed on to \code{\link[metafor]{rma.uni}} or
#'   \code{\link[metafor]{rma.glmm}}.
#' @param \dots Additional arguments passed on to
#'   \code{\link[metafor]{rma.glmm}} function and to catch deprecated
#'   arguments.
#' 
#' @details
#' Calculation of common and random effects estimates for meta-analyses
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
#' \item Vaccine efficacy or vaccine effectiveness (\code{sm = "VE"})
#' }
#'
#' Note, mathematically, odds ratios and diagnostic odds ratios are
#' identical, however, the labels in printouts and figures
#' differ. Furthermore, log risk ratio (logRR) and log vaccine ratio
#' (logVR) are mathematical identical, however, back-transformed
#' results differ as vaccine efficacy or effectiveness is defined as
#' \code{VE = 100 * (1 - RR)}.
#' 
#' A three-level random effects meta-analysis model (Van den Noortgate
#' et al., 2013) is utilised if argument \code{cluster} is used and at
#' least one cluster provides more than one estimate. Internally,
#' \code{\link[metafor]{rma.mv}} is called to conduct the analysis and
#' \code{\link[metafor]{weights.rma.mv}} with argument \code{type =
#' "rowsum"} is used to calculate random effects weights.
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
#' By default, both common effect (also called common effect) and
#' random effects models are considered (see arguments \code{common}
#' and \code{random}). If \code{method} is \code{"MH"} (default), the
#' Mantel-Haenszel method (Greenland & Robins, 1985; Robins et al.,
#' 1986) is used to calculate the common effect estimate; if
#' \code{method} is \code{"Inverse"}, inverse variance weighting is
#' used for pooling (Fleiss, 1993); if \code{method} is \code{"Peto"},
#' the Peto method is used for pooling (Yusuf et al., 1985); if
#' \code{method} is \code{"SSW"}, the sample size method is used for
#' pooling (Bakbergenuly et al., 2020).
#'
#' While the Mantel-Haenszel and Peto method are defined under the
#' common effect model, random effects variants based on these methods
#' are also implemented in \code{metabin}. Following RevMan 5, the
#' Mantel-Haenszel estimator is used in the calculation of the
#' between-study heterogeneity statistic Q which is used in the
#' DerSimonian-Laird estimator (DerSimonian and Laird,
#' 1986). Accordingly, the results for the random effects
#' meta-analysis using the Mantel-Haenszel or inverse variance method
#' are typically very similar. For the Peto method, Peto's log odds
#' ratio, i.e. \code{(O-E) / V} and its standard error \code{sqrt(1 /
#' V)} with \code{O-E} and \code{V} denoting "Observed minus Expected"
#' and its variance, are utilised in the random effects
#' model. Accordingly, results of a random effects model using
#' \code{sm = "Peto"} can be different to results from a random
#' effects model using \code{sm = "MH"} or \code{sm =
#' "Inverse"}. Note, the random effects estimate is based on the
#' inverse variance method for all methods discussed so far.
#' 
#' A distinctive and frequently overlooked advantage of binary
#' endpoints is that individual patient data (IPD) can be extracted
#' from a two-by-two table.  Accordingly, statistical methods for IPD,
#' i.e., logistic regression and generalised linear mixed models, can
#' be utilised in a meta-analysis of binary outcomes (Stijnen et al.,
#' 2010; Simmonds et al., 2016).
#' 
#' R package \bold{brglm2} must be available to fit a one-stage logistic
#' regression model with penalised likelihood (Evrenoglou et al., 2022).
#' The estimation of the summary odds ratio relies on the maximisation of the 
#' likelihood function, penalised using a Firth-type correction. This
#' penalisation aims to reduce bias in cases with rare events and a small
#' number of available studies. However, this method is not restricted 
#' to only such cases and can be applied more generally to binary data. Note,
#' with this type of penalisation, all studies can be included in the analysis,
#' regardless of the total number of observed events. This allows both single
#' and double zero studies to be included without any continuity correction.
#' The random effects model uses a multiplicative heterogeneity parameter
#' \eqn{\phi}, added to the model as an \emph{ad hoc} term. The estimation of
#' this parameter relies on a modified expression of Pearson's statistic, which
#' accounts for sparse data. An estimate of \eqn{\phi} equal to 1 indicates the
#' absence of heterogeneity.
#' 
#' Generalised linear mixed models are available
#' (argument \code{method = "GLMM"}) for the odds ratio as summary measure for
#' the common effect and random effects model by calling the
#' \code{\link[metafor]{rma.glmm}} function from R package
#' \bold{metafor} internally. 
#'
#' Four different GLMMs are available for
#' meta-analysis with binary outcomes using argument \code{model.glmm}
#' (which corresponds to argument \code{model} in the
#' \code{\link[metafor]{rma.glmm}} function):
#' \tabular{cl}{
#' 1. \tab Logistic regression model with common study effects
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
#' can be provided using argument '\code{\dots}' in \code{metabin} are
#' described in \code{\link[metafor]{rma.glmm}} where you can also
#' find information on the iterative algorithms used for estimation.
#' Note, regardless of which value is used for argument
#' \code{model.glmm}, results for two different GLMMs are calculated:
#' common effect model (with fixed treatment effect) and random
#' effects model (with random treatment effects).
#' }
#' 
#' \subsection{Continuity correction}{
#'
#' Four approaches are available to apply a continuity correction:
#' \itemize{
#' \item Only studies with a zero cell count (\code{method.incr =
#'   "only0"})
#' \item All studies if at least one study has a zero cell count
#'   (\code{method.incr = "if0all"})
#' \item All studies irrespective of zero cell counts
#'   (\code{method.incr = "all"})
#' \item Use values provided in arguments \code{incr.e} and \code{incr.c}
#'   (\code{method.incr = "user"})
#' }
#'
#' By default, a continuity correction is only applied to studies with
#' a zero cell count (\code{method.incr = "only0"}). This method
#' showed the best performance for the odds ratio in a simulation
#' study under the random effects model (Weber et al., 2020).
#'
#' The continuity correction method is used both to calculate
#' individual study results with confidence limits and to conduct
#' meta-analysis based on the inverse variance method. For the risk
#' difference, the method is only considered to calculate standard
#' errors and confidence limits.  For Peto method and GLMMs no
#' continuity correction is used in the meta-analysis. Furthermore,
#' the continuity correction is ignored for individual studies for the
#' Peto method.
#'
#' For studies with a zero cell count, by default, 0.5 (argument
#' \code{incr}) is added to all cell frequencies for the odds ratio or
#' only the number of events for the risk ratio (argument
#' \code{RR.Cochrane = FALSE}, default). The increment is added to all
#' cell frequencies for the risk ratio if argument \code{RR.Cochrane =
#' TRUE}. For the risk difference, \code{incr} is only added to all
#' cell frequencies to calculate the standard error. Finally, a
#' treatment arm continuity correction is used if \code{incr = "TACC"}
#' (Sweeting et al., 2004; Diamond et al., 2007).
#'
#' For odds ratio and risk ratio, treatment estimates and standard
#' errors are only calculated for studies with zero or all events in
#' both groups if \code{allstudies = TRUE}.
#'
#' For the Mantel-Haenszel method, by default (if \code{MH.exact} is
#' FALSE), \code{incr} is added to cell frequencies of a study with a
#' zero cell count in the calculation of the pooled risk ratio or odds
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
#' \subsection{Subgroup analysis}{
#' 
#' Argument \code{subgroup} can be used to conduct subgroup analysis for
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
#' Internally, both common effect and random effects models are
#' calculated regardless of values choosen for arguments
#' \code{common} and \code{random}. Accordingly, the estimate
#' for the random effects model can be extracted from component
#' \code{TE.random} of an object of class \code{"meta"} even if
#' argument \code{random = FALSE}. However, all functions in R
#' package \bold{meta} will adequately consider the values for
#' \code{common} and \code{random}. E.g. function
#' \code{\link{print.meta}} will not print results for the random
#' effects model if \code{random = FALSE}.
#'
#' A prediction interval will only be shown if \code{prediction =
#' TRUE}.
#' }
#' 
#' @return
#' An object of class \code{c("metabin", "meta")} with corresponding
#' generic functions (see \code{\link{meta-object}}).
#'
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{update.meta}},
#'   \code{\link{forest}}, \code{\link{funnel}},
#'   \code{\link{metabias}}, \code{\link{metacont}},
#'   \code{\link{metagen}}, \code{\link{metareg}},
#'   \code{\link{print.meta}}
#' 
#' @references
#' Bakbergenuly I, Hoaglin DC, Kulinskaya E (2020):
#' Methods for estimating between-study variance and overall
#' effect in meta-analysis of odds-ratios.
#' \emph{Research Synthesis Methods},
#' \bold{11}, 426--42
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
#' Evrenoglou T, White IR, Afach S, Mavridis D, Chaimani A. (2022):
#' Network meta-analysis of rare events using penalized likelihood regression.
#' \emph{Statistics in Medicine},
#' \bold{41}, 5203--19
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
#' \emph{Review Manager (RevMan)} [Computer program]. Version 5.4.
#' The Cochrane Collaboration, 2020
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
#' Van den Noortgate W, López-López JA, Marín-Martínez F, Sánchez-Meca J (2013):
#' Three-level meta-analysis of dependent effect sizes.
#' \emph{Behavior Research Methods},
#' \bold{45}, 576--94
#' 
#' Viechtbauer W (2010):
#' Conducting meta-analyses in R with the metafor package.
#' \emph{Journal of Statistical Software},
#' \bold{36}, 1--48
#' 
#' Weber F, Knapp G, Ickstadt K, Kundt G, Glass Ä (2020):
#' Zero-cell corrections in random-effects meta-analyses.
#' \emph{Research Synthesis Methods},
#' \bold{11}, 913--9
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
#'   data = Olkin1995, subset = c(41, 47, 51, 59),
#'   studlab = paste(author, year),
#'   method = "Inverse")
#' m1
#' # Show results for individual studies
#' summary(m1)
#' 
#' # Use different subset of Olkin (1995)
#' #
#' m2 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'   data = Olkin1995, subset = year < 1970,
#'   studlab = paste(author, year),
#'   method = "Inverse")
#' m2
#' forest(m2)
#' 
#' # Meta-analysis with odds ratio as summary measure
#' #
#' m3 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'   data = Olkin1995, subset = year < 1970,
#'   studlab = paste(author, year),
#'   sm = "OR", method = "Inverse")
#' # Same meta-analysis result using 'update.meta' function
#' m3 <- update(m2, sm = "OR")
#' m3
#' 
#' # Meta-analysis based on Mantel-Haenszel method (with odds ratio as
#' # summary measure)
#' #
#' m4 <- update(m3, method = "MH")
#' m4
#' 
#' # Meta-analysis based on Peto method (only available for odds ratio
#' # as summary measure)
#' #
#' m5 <- update(m3, method = "Peto")
#' m5
#' 
#' \dontrun{
#' # Meta-analyses using generalised linear mixed models (GLMM)
#' 
#' # Logistic regression model with (k = 4) fixed study effects
#' # (default: model.glmm = "UM.FS")
#' m6 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'   studlab = paste(author, year),
#'   data = Olkin1995, subset = year < 1970, method = "GLMM")
#' # Same results:
#' m6 <- update(m2, method = "GLMM")
#' m6
#' 
#' # Mixed-effects logistic regression model with random study effects
#' m7 <- update(m6, model.glmm = "UM.RS")
#' #
#' # Use additional argument 'nAGQ' for internal call of 'rma.glmm'
#' # function
#' #
#' m7 <- update(m6, model.glmm = "UM.RS", nAGQ = 1)
#' m7
#' 
#' # Generalised linear mixed model (conditional Hypergeometric-Normal)
#' # (R package 'BiasedUrn' must be available)
#' if (requireNamespace("BiasedUrn", quietly = TRUE)) {
#'  m8 <- update(m6, model.glmm = "CM.EL")
#'  m8
#' }
#' 
#' # Generalised linear mixed model (conditional Binomial-Normal)
#' m9 <- update(m6, model.glmm = "CM.AL")
#' m9
#' 
#' # Logistic regression model with (k = 70) fixed study effects
#' m10 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
#'    studlab = paste(author, year),
#'    data = Olkin1995, method = "GLMM")
#' m10
#' 
#' # Mixed-effects logistic regression model with random study effects
#' update(m10, model.glmm = "UM.RS")
#' 
#' # Conditional Hypergeometric-Normal GLMM (with long computation time)
#' system.time(m11 <- update(m10, model.glmm = "CM.EL"))
#' m11
#' 
#' # Generalised linear mixed model (conditional Binomial-Normal)
#' update(m10, model.glmm = "CM.AL")
#' }
#' 
#' @export metabin


metabin <- function(event.e, n.e, event.c, n.c, studlab,
                    ##
                    data = NULL, subset = NULL, exclude = NULL,
                    cluster = NULL, rho = 0,
                    #
                    weights = NULL,
                    weights.common = weights, weights.random = weights,
                    #
                    method = ifelse(tau.common, "Inverse", gs("method")),
                    sm =
                      ifelse(!is.na(charmatch(tolower(method),
                                              c("peto", "glmm", "lrp", "ssw"),
                                              nomatch = NA)),
                             "OR", gs("smbin")),
                    incr = gs("incr"), method.incr = gs("method.incr"),
                    allstudies = gs("allstudies"),
                    incr.e = if (length(incr) > 1) incr else NULL,
                    incr.c = if (length(incr) > 1) incr else NULL,
                    #
                    level = gs("level"),
                    ##
                    MH.exact = gs("MH.exact"), RR.Cochrane = gs("RR.Cochrane"),
                    Q.Cochrane =
                      gs("Q.Cochrane") & method == "MH" & method.tau == "DL",
                    model.glmm = gs("model.glmm"),
                    ##
                    common = gs("common"),
                    random = gs("random") | !is.null(tau.preset),
                    overall = common | random,
                    overall.hetstat =
                      if (is.null(gs("overall.hetstat")))
                        common | random
                      else
                        gs("overall.hetstat"),   
                    prediction = gs("prediction") | !missing(method.predict),
                    ##
                    method.tau,
                    method.tau.ci = gs("method.tau.ci"),
                    level.hetstat = gs("level.hetstat"),
                    tau.preset = NULL, TE.tau = NULL,
                    tau.common = gs("tau.common"),
                    detail.tau = NULL,
                    #
                    method.I2 = gs("method.I2"),
                    #
                    level.ma = gs("level.ma"),
                    method.common.ci = gs("method.common.ci"),
                    method.random.ci = gs("method.random.ci"),
                    adhoc.hakn.ci = gs("adhoc.hakn.ci"),
                    ##
                    level.predict = gs("level.predict"),
                    method.predict = gs("method.predict"),
                    adhoc.hakn.pi = gs("adhoc.hakn.pi"),
                    seed.predict = NULL,
                    ##
                    method.bias = ifelse(sm == "OR", "Harbord",
                                  ifelse(sm == "DOR", "Deeks",
                                         gs("method.bias"))),
                    ##
                    backtransf = gs("backtransf"),
                    pscale = 1,
                    ##
                    text.common = gs("text.common"),
                    text.random = gs("text.random"),
                    text.predict = gs("text.predict"),
                    text.w.common = gs("text.w.common"),
                    text.w.random = gs("text.w.random"),
                    ##
                    title = gs("title"), complab = gs("complab"),
                    outclab = "",
                    #
                    label.e = gs("label.e"), label.c = gs("label.c"),
                    label.left = gs("label.left"),
                    label.right = gs("label.right"),
                    col.label.left = gs("col.label.left"),
                    col.label.right = gs("col.label.right"),
                    #
                    subgroup, subgroup.name = NULL,
                    print.subgroup.name = gs("print.subgroup.name"),
                    sep.subgroup = gs("sep.subgroup"),
                    test.subgroup = gs("test.subgroup"),
                    prediction.subgroup = gs("prediction.subgroup"),
                    seed.predict.subgroup = NULL,
                    ##
                    byvar, hakn, adhoc.hakn,
                    ##
                    print.CMH = gs("print.CMH"),
                    ##
                    keepdata = gs("keepdata"),
                    warn = gs("warn"), warn.deprecated = gs("warn.deprecated"),
                    ##
                    control = NULL,
                    ...
                    ) {


  ##
  ##
  ## (1) Check arguments
  ##
  ##
  
  args <- list(...)
  nam.args <- names(args)
  #
  missing.sm <- missing(sm)
  missing.subgroup <- missing(subgroup)
  missing.byvar <- missing(byvar)
  missing.overall <- missing(overall)
  missing.overall.hetstat <- missing(overall.hetstat)
  missing.test.subgroup <- missing(test.subgroup)
  #
  missing.event.c <- missing(event.c)
  missing.n.e <- missing(n.e)
  missing.n.c <- missing(n.c)
  #
  missing.studlab <- missing(studlab)
  #
  missing.incr <- missing(incr)
  avail.incr <- !missing.incr && !is.null(incr)
  #
  missing.method.incr <- missing(method.incr)
  avail.method.incr <- !missing.method.incr && !is.null(method.incr)
  #
  missing.allstudies <- missing(allstudies)
  #
  missing.incr.e <- missing(incr.e)
  missing.incr.c <- missing(incr.c)
  #
  missing.method.tau <- missing(method.tau)
  missing.tau.common <- missing(tau.common)
  missing.method.predict <- missing(method.predict)
  missing.method <- missing(method)
  missing.Q.Cochrane <- missing(Q.Cochrane)
  missing.level.ma <- missing(level.ma)
  missing.common <- missing(common)
  missing.random <- missing(random)
  missing.method.common.ci <- missing(method.common.ci)
  missing.method.random.ci <- missing(method.random.ci)
  missing.method.I2 <- missing(method.I2)
  #
  missing.hakn <- missing(hakn)
  missing.adhoc.hakn.ci <- missing(adhoc.hakn.ci)
  missing.adhoc.hakn <- missing(adhoc.hakn)
  missing.RR.Cochrane <- missing(RR.Cochrane)
  #
  missing.subgroup.name <- missing(subgroup.name)
  missing.print.subgroup.name <- missing(print.subgroup.name)
  missing.sep.subgroup <- missing(sep.subgroup)
  #
  missing.label.e <- missing(label.e)
  missing.label.c <- missing(label.c)
  missing.complab <- missing(complab)
  #
  missing.cluster <- missing(cluster)
  #
  chknumeric(rho, min = -1, max = 1)
  ##
  chknull(sm)
  sm <- setchar(sm, gs("sm4bin"))
  #
  chklevel(level)
  #
  method <- setchar(method, gs("meth4bin"))
  #
  is.glmm <- method == "GLMM"
  is.lrp <- method == "LRP"
  #
  method.common.ci <- setchar(method.common.ci, gs("meth4common.ci"))
  #
  if (method != "Inverse" & method.common.ci == "IVhet") {
    if (!missing.method.common.ci)
      warning("Argument 'method.common.ci = \"IVhet\"' only available ",
              "if 'method = \"Inverse\".",
              call. = FALSE)
    method.common.ci <- "classic"
  }
  #
  if (missing.method.tau) {
    if (is.lrp)
      method.tau <- "DL"
    else if (is.glmm)
      method.tau <- "ML"
    else
      method.tau <- gs("method.tau")
  }
  #
  method.tau <- setchar(method.tau, c(gs("meth4tau"), "KD"))
  ##
  tau.common <- replaceNULL(tau.common, FALSE)
  chklogical(tau.common)
  #
  method.I2 <- setchar(method.I2, gs("meth4i2"))
  #
  chklogical(prediction)
  chklevel(level.predict)
  ##
  method.predict <- setchar(method.predict, gs("meth4pi"))
  ##
  method.tau <-
    set_method_tau(method.tau, missing.method.tau,
                   method.predict, missing.method.predict)
  method.predict <-
    set_method_predict(method.predict, missing.method.predict,
                       method.tau, missing.method.tau)
  ##
  if (any(method.predict == "NNF"))
    is_installed_package("pimeta", argument = "method.predict", value = "NNF")
  ##
  adhoc.hakn.pi <- setchar(replaceNA(adhoc.hakn.pi, ""), gs("adhoc4hakn.pi"))
  #
  chklogical(backtransf)
  ##
  chknumeric(pscale, length = 1)
  ##
  if (!is.null(text.common))
    chkchar(text.common, length = 1)
  if (!is.null(text.random))
    chkchar(text.random)
  if (!is.null(text.predict))
    chkchar(text.predict)
  if (!is.null(text.w.common))
    chkchar(text.w.common, length = 1)
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
      warning("Argument 'pscale' only considered for risk differences.",
              call. = FALSE)
    pscale <- 1
  }
  #
  method.incr <- setchar(method.incr, gs("meth4incr"))
  ##
  chklogical(allstudies)
  #
  chklogical(MH.exact)
  chklogical(Q.Cochrane)
  if (Q.Cochrane & (method != "MH" | method.tau != "DL")) {
    warning("Argument 'Q.Cochrane' only considered for ",
            "Mantel-Haenszel method in combination with ",
            "DerSimonian-Laird estimator.",
            call. = FALSE)
    Q.Cochrane <- FALSE
  }
  ##
  if (length(model.glmm) == 0)
    model.glmm <- gs("model.glmm")
  model.glmm <- setchar(model.glmm, c("UM.FS", "UM.RS", "CM.EL", "CM.AL", ""))
  #
  chklogical(print.CMH)
  ##
  if (sm == "ASD") {
    method <- "Inverse"
    if (!missing.Q.Cochrane && Q.Cochrane)
      warning("Argument 'Q.Cochrane' only considered for ",
              "Mantel-Haenszel method in combination with ",
              "DerSimonian-Laird estimator.",
              call. = FALSE)
    Q.Cochrane <- FALSE
  }
  ##
  ## Check for deprecated arguments in '...'
  ##
  chklogical(warn.deprecated)
  ##
  level.ma <- deprecated(level.ma, missing.level.ma, args, "level.comb",
                         warn.deprecated)
  chklevel(level.ma)
  ##
  common <- deprecated(common, missing.common, args, "comb.fixed",
                       warn.deprecated)
  common <- deprecated(common, missing.common, args, "fixed",
                       warn.deprecated)
  chklogical(common)
  ##
  random <- deprecated(random, missing.random, args, "comb.random",
                       warn.deprecated)
  chklogical(random)
  ##
  method.random.ci <-
    deprecated2(method.random.ci, missing.method.random.ci,
                hakn, missing.hakn,
                warn.deprecated)
  if (is.logical(method.random.ci))
    if (method.random.ci)
      method.random.ci <- "HK"
    else
      method.random.ci <- "classic"
  method.random.ci <- setchar(method.random.ci, gs("meth4random.ci"))
  ##
  adhoc.hakn.ci <-
    deprecated2(adhoc.hakn.ci, missing.adhoc.hakn.ci,
                adhoc.hakn, missing.adhoc.hakn, warn.deprecated)
  adhoc.hakn.ci <- setchar(replaceNA(adhoc.hakn.ci, ""), gs("adhoc4hakn.ci"))
  #
  missing.subgroup.name <- missing.subgroup.name
  subgroup.name <-
    deprecated(subgroup.name, missing.subgroup.name, args, "bylab",
               warn.deprecated)
  ##
  print.subgroup.name <-
    deprecated(print.subgroup.name, missing.print.subgroup.name,
               args, "print.byvar", warn.deprecated)
  print.subgroup.name <-
    replaceNULL(print.subgroup.name, gs("print.subgroup.name"))
  chklogical(print.subgroup.name)
  ##
  sep.subgroup <-
    deprecated(sep.subgroup, missing.sep.subgroup, args, "byseparator",
               warn.deprecated)
  if (!is.null(sep.subgroup))
    chkchar(sep.subgroup, length = 1)
  ##
  RR.Cochrane <-
    deprecated(RR.Cochrane, missing.RR.Cochrane, args, "RR.cochrane",
               warn.deprecated)
  chklogical(RR.Cochrane)
  ##
  ## Some more checks
  ##
  chklogical(overall)
  chklogical(overall.hetstat)
  #
  # Ignore deprecated arguments 'addincr' and 'allincr'
  #
  txt.ignore <- "(deprecated); use argument 'method.incr'"
  #
  if (!is.na(charmatch("addincr", nam.args)))
    warn_ignore_input(addincr, TRUE, txt.ignore)
  if (!is.na(charmatch("allincr", nam.args)))
    warn_ignore_input(allincr, TRUE, txt.ignore)
  
  
  ##
  ##
  ## (2) Read data
  ##
  ##
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (nulldata) {
    data <- sfsp
    data.pairwise <- FALSE
  }
  else
    data.pairwise <- inherits(data, "pairwise")
  #
  # Catch 'event.e', 'n.e', 'event.c', 'n.c', 'studlab', 'subgroup', 'incr',
  # 'incr.e' and 'incr.c' from data:
  #
  event.e <- catch("event.e", mc, data, sfsp)
  chknull(event.e)
  #
  if (inherits(event.e, "pairwise")) {
    is.pairwise <- TRUE
    #
    type <- attr(event.e, "type")
    if (type != "binary")
      stop("Wrong type for pairwise() object: '", type, "'.", call. = FALSE)
    #
    txt.ignore <- "as first argument is a pairwise object"
    #
    warn_ignore_input(event.c, !missing.event.c, txt.ignore)
    warn_ignore_input(n.e, !missing.n.e, txt.ignore)
    warn_ignore_input(n.c, !missing.n.c, txt.ignore)
    warn_ignore_input(incr.e, !missing.incr.e, txt.ignore)
    warn_ignore_input(incr.c, !missing.incr.c, txt.ignore)
    #
    warn_ignore_input(studlab, !missing.studlab, txt.ignore)
    #
    if (missing.sm)
      sm <- attr(event.e, "sm")
    #
    if (missing.method)
      method <- attr(event.e, "method")
    #
    if (!avail.method.incr & !avail.incr)
      method.incr <- "user"
    #
    missing.sm <- FALSE
    missing.method <- FALSE
    #
    avail.method.incr <- TRUE
    missing.incr <- FALSE
    missing.method.incr <- FALSE
    #
    reference.group <- attr(event.e, "reference.group")
    #
    studlab <- event.e$studlab
    #
    treat1 <- event.e$treat1
    treat2 <- event.e$treat2
    #
    event.c <- event.e$event2
    #
    n.e <- event.e$n1
    n.c <- event.e$n2
    #
    incr.e <- event.e$incr1
    incr.c <- event.e$incr2
    #
    if (avail.incr | method.incr != "user") {
      incr.e <- NULL
      incr.c <- NULL
    }
    #
    if (!avail.incr)
      incr <- attr(event.e, "incr")
    #
    pairdata <- event.e
    data <- event.e
    nulldata <- FALSE
    #
    event.e <- event.e$event1
    #
    wo <- treat1 == reference.group
    #
    if (any(wo)) {
      ttreat1 <- treat1
      treat1[wo] <- treat2[wo]
      treat2[wo] <- ttreat1[wo]
      #
      tevent.e <- event.e
      event.e[wo] <- event.c[wo]
      event.c[wo] <- tevent.e[wo]
      #
      tn.e <- n.e
      n.e[wo] <- n.c[wo]
      n.c[wo] <- tn.e[wo]
      #
      if (!(is.null(incr.e) | is.null(incr.c))) {
        tincr.e <- incr.e
        incr.e[wo] <- incr.c[wo]
        incr.c[wo] <- tincr.e[wo]
      }
    }
    #
    if (missing.subgroup & missing.byvar) {
      subgroup <- paste(treat1, treat2, sep = " vs ")
      #
      if (length(unique(subgroup)) == 1) {
        if (missing.complab)
          complab <- unique(subgroup)
        #
        if (missing.label.e)
          label.e <- unique(treat1)
        if (missing.label.c)
          label.c <- unique(treat2)
        #
        subgroup <- NULL
      }
      else {
        if (missing.overall)
          overall <- FALSE
        if (missing.overall.hetstat)
          overall.hetstat <- FALSE
        if (missing.test.subgroup)
          test.subgroup <- FALSE
      }
    }
    else {
      subgroup <- catch("subgroup", mc, data, sfsp)
      byvar <- catch("byvar", mc, data, sfsp)
      #
      subgroup <- deprecated2(subgroup, missing.subgroup, byvar, missing.byvar,
                              warn.deprecated)
    }
  }
  else {
    is.pairwise <- FALSE
    #
    if (missing.sm && !is.null(data) && !is.null(attr(data, "sm")))
      sm <- attr(data, "sm")
    #
    n.e <- catch("n.e", mc, data, sfsp)
    #
    event.c <- catch("event.c", mc, data, sfsp)
    n.c <- catch("n.c", mc, data, sfsp)
    #
    studlab <- catch("studlab", mc, data, sfsp)
    #
    subgroup <- catch("subgroup", mc, data, sfsp)
    byvar <- catch("byvar", mc, data, sfsp)
    #
    subgroup <- deprecated2(subgroup, missing.subgroup, byvar, missing.byvar,
                            warn.deprecated)
    #
    if (!missing.incr)
      incr <- catch("incr", mc, data, sfsp)
    #
    if (!missing.incr.e)
      incr.e <- catch("incr.e", mc, data, sfsp)
    if (!missing.incr.c)
      incr.c <- catch("incr.c", mc, data, sfsp)
  }
  #
  method.bias <- setmethodbias(method.bias)
  #
  is.tacc <- FALSE
  #
  if (is.numeric(incr))
    chknumeric(incr, min = 0)
  else {
    incr <- setchar(incr, "TACC",
                    "should be numeric or the character string \"TACC\"")
    is.tacc <- TRUE
  }
  #
  avail.incr.e <- !is.null(incr.e)
  avail.incr.c <- !is.null(incr.c)
  avail.incr.both <- avail.incr.e & avail.incr.c
  #
  if (avail.incr.e + avail.incr.c == 1)
    stop("Arguments 'incr.e' and 'incr.c' are required together.",
         call. = FALSE)
  #
  if (avail.incr.e)
    chknumeric(incr.e, min = 0)
  #
  if (avail.incr.c)
    chknumeric(incr.c, min = 0)
  #
  if (avail.incr.both) {
    if (!avail.method.incr)
      method.incr <- "user"
    #
    txt.ignore <- "as arguments 'incr.e' and 'incr.c' are provided"
    #
    warn_set_input(method.incr, method.incr != "user", '"user"', txt.ignore)
    #
    if (length(incr) != length(incr.e) | any(incr != incr.e))
      warn_ignore_input(incr, avail.incr, txt.ignore)
  }
  else {
    addincr <- allincr <- FALSE
    #
    if (method.incr == "all")
      addincr <- TRUE
    else if (method.incr == "if0all")
      allincr <- TRUE
    addincr <- allincr <- FALSE
    #
    if (!(sm == "ASD" | method %in% c("Peto", "GLMM"))) {
      if (method.incr == "all")
        addincr <- TRUE
      else if (method.incr == "if0all")
        allincr <- TRUE
    }
  }
  #
  k.All <- length(event.e)
  #
  chknull(n.e)
  chknull(event.c)
  chknull(n.c)
  #
  studlab <- setstudlab(studlab, k.All)
  #
  by <- !is.null(subgroup)
  #
  # Catch 'subset', 'exclude' and 'cluster' from data:
  #
  subset <- catch("subset", mc, data, sfsp)
  missing.subset <- is.null(subset)
  ##
  exclude <- catch("exclude", mc, data, sfsp)
  missing.exclude <- is.null(exclude)
  ##
  cluster <- catch("cluster", mc, data, sfsp)
  with.cluster <- !is.null(cluster)
  #
  # Catch 'weights', 'weights.common', and 'weights.random' from data:
  #
  if (!missing(weights))
    weights <- catch("weights", mc, data, sfsp)
  if (!missing(weights.common))
    weights.common <- catch("weights.common", mc, data, sfsp)
  if (!missing(weights.random))
    weights.random <- catch("weights.random", mc, data, sfsp)
  #
  if (!is.null(weights) & is.null(weights.common))
    weights.common <- weights
  #
  if (!is.null(weights) & is.null(weights.random))
    weights.random <- weights
  #
  usw.common <- !is.null(weights.common)
  usw.random <- !is.null(weights.random)
  #
  if (usw.common)
    chknumeric(weights.common, min = 0)
  #
  if (usw.random)
    chknumeric(weights.random, min = 0)
  #
  if (usw.common & method != "Inverse")
    stop("User-specified weights for the common effect model only implemented ",
         "for the inverse variance method (method = \"Inverse\").",
         call. = FALSE)
  #
  if (usw.random & method %in% c("GLMM", "LRP", "SSW"))
    stop("User-specified weights for the random effects model not implemented ",
         "for method = \"GLMM\",  \"LRP\" or  \"SSW\".",
         call. = FALSE)
  #
  # Check variable values
  #
  chknumeric(event.e)
  chknumeric(n.e)
  chknumeric(event.c)
  chknumeric(n.c)
  #
  # Recode integer as numeric:
  #
  event.e <- int2num(event.e)
  n.e     <- int2num(n.e)
  event.c <- int2num(event.c)
  n.c     <- int2num(n.c)
  
  
  ##
  ##
  ## (3) Check length of essential variables
  ##
  ##
  
  chklength(n.e, k.All, fun)
  chklength(event.c, k.All, fun)
  chklength(n.c, k.All, fun)
  chklength(studlab, k.All, fun)
  #
  if (with.cluster)
    chklength(cluster, k.All, fun)
  #
  if (usw.common) {
    if (length(weights.common) == 1)
      weights.common <- rep(weights.common, k.All)
    else
      chklength(weights.common, k.All, fun)
  }
  #
  if (usw.random) {
    if (length(weights.random) == 1)
      weights.random <- rep(weights.random, k.All)
    else
      chklength(weights.random, k.All, fun)
  }
  #
  if (length(incr) > 1)
    chklength(incr, k.All, fun)
  #
  if (avail.incr.e) {
    if (length(incr.e) == 1)  
      incr.e <- rep_len(incr.e, k.All)
    else
      chklength(incr.e, k.All, fun)
  }
  #
  if (avail.incr.c) {
    if (length(incr.c) == 1)  
      incr.c <- rep_len(incr.c, k.All)
    else
      chklength(incr.c, k.All, fun)
  }
  #
  if (by) {
    chklength(subgroup, k.All, fun)
    chklogical(test.subgroup)
    chklogical(prediction.subgroup)
  }
  ##
  ## Additional checks
  ##
  if (!by & tau.common) {
    if (warn)
      warning("Value for argument 'tau.common' set to FALSE as ",
              "argument 'subgroup' is missing.",
              call. = FALSE)
    tau.common <- FALSE
  }
  #
  if (by & !tau.common & !is.null(tau.preset)) {
    if (warn)
      warning("Argument 'tau.common' set to TRUE as ",
              "argument tau.preset is not NULL.",
              call. = FALSE)
    tau.common <- TRUE
  }
  #
  if (method.incr == "user" & !avail.incr.both) {
    if (is.null(incr.e) | is.null(incr.c))
      stop("Non-null input for arguments 'incr.e' and 'incr.c' required if ",
           "'method.incr = \"user\"'.",
           call. = FALSE)
    else
      stop("Arguments 'incr.e' and 'incr.c' are required if ",
           "'method.incr = \"user\"'.",
           call. = FALSE)
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
      data <- data.frame(.studlab = studlab)
    else
      data$.studlab <- studlab
    ##
    data$.event.e <- event.e
    data$.n.e <- n.e
    data$.event.c <- event.c
    data$.n.c <- n.c
    #
    data$.incr.e <- NA
    data$.incr.c <- NA
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
    ##
    if (with.cluster)
      data$.id <- data$.cluster <- cluster
    #
    if (usw.common)
      data$.weights.common <- weights.common
    #
    if (usw.random)
      data$.weights.random <- weights.random
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
    cluster <- cluster[subset]
    exclude <- exclude[subset]
    #
    weights.common <- weights.common[subset]
    weights.random <- weights.random[subset]
    #
    incr.e <- incr.e[subset]
    incr.c <- incr.c[subset]
    #
    if (length(incr) > 1)
      incr <- incr[subset]
    #
    if (by)
      subgroup <- subgroup[subset]
  }
  #
  if (missing.subgroup & is.pairwise & by) {
    if (length(unique(subgroup)) == 1) {
      by <- FALSE
      #
      if (missing.complab)
        complab <- unique(subgroup)
      #
      subgroup <- NULL
      #
      if (keepdata)
        data$.subgroup <- NULL
      #
      if (missing.overall)
        overall <- TRUE
      if (missing.overall.hetstat)
        overall.hetstat <- TRUE
    }
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
    common <- FALSE
    random <- FALSE
    prediction <- FALSE
    overall <- FALSE
    overall.hetstat <- FALSE
  }
  ##
  if (by) {
    chkmiss(subgroup)
    ##
    if (missing.subgroup.name & is.null(subgroup.name)) {
      if (!missing.subgroup)
        subgroup.name <- byvarname("subgroup", mc)
      else if (!missing.byvar)
        subgroup.name <- byvarname("byvar", mc)
    }
  }
  ##
  if (!is.null(subgroup.name))
    chkchar(subgroup.name, length = 1)
  
  
  #
  #
  # (7) Continuity correction
  #
  #
  
  #
  # Include non-informative studies?
  # (i.e. studies with either zero or all events in both groups)
  #
  allevents <- event.c == n.c & event.e == n.e
  #
  if (sm == "RD" | sm == "ASD")
    incl <- rep(1, k.all)
  else {
    if (allstudies)
      incl <- rep(1, k.all)
    else {
      if (sm %in% c("OR", "DOR"))
        incl <- ifelse((event.c == 0   & event.e == 0) |
                         (event.c == n.c & event.e == n.e), NA, 1)
      if (sm %in% c("RR", "VE"))
        incl <- ifelse((event.c == 0 & event.e == 0), NA, 1)
    }
  }
  #
  # Exclude studies from meta-analysis:
  #
  sel1 <- event.e > n.e
  sel2 <- event.c > n.c
  if ((any(sel1, na.rm = TRUE)) & warn)
    warning("Studies with event.e > n.e get no weight in meta-analysis.",
            call. = FALSE)
  if ((any(sel2, na.rm = TRUE)) & warn)
    warning("Studies with event.c > n.c get no weight in meta-analysis.",
            call. = FALSE)
  incl[sel1 | sel2] <- NA
  #
  sel3 <- n.e <= 0 | n.c <= 0
  if ((any(sel3, na.rm = TRUE)) & warn)
    warning("Studies with non-positive values for n.e and / or n.c ",
            "get no weight in meta-analysis.",
            call. = FALSE)
  incl[sel3] <- NA
  #
  sel4 <- event.e < 0 | event.c < 0
  if ((any(sel4, na.rm = TRUE)) & warn)
    warning("Studies with negative values for event.e and / or event.c ",
            "get no weight in meta-analysis.",
            call. = FALSE)
  incl[sel4] <- NA
  #
  # Sparse computation
  #
  sel <- switch(sm,
                OR = ((n.e - event.e) == 0 | event.e == 0 |
                        (n.c - event.c) == 0 | event.c == 0),
                RD = ((n.e - event.e) == 0 | event.e == 0 |
                        (n.c - event.c) == 0 | event.c == 0),
                RR = ((n.e - event.e) == 0 | event.e == 0 |
                        (n.c - event.c) == 0 | event.c == 0),
                VE = ((n.e - event.e) == 0 | event.e == 0 |
                        (n.c - event.c) == 0 | event.c == 0),
                ASD = rep(FALSE, length(event.e)),
                DOR = ((n.e - event.e) == 0 | event.e == 0 |
                         (n.c - event.c) == 0 | event.c == 0))
  #
  sel[is.na(incl)] <- FALSE
  incl[is.na(incl)] <- 0
  #
  sparse <- any(sel, na.rm = TRUE)
  #
  # Check for pairwise comparisons with zero cell frequencies in both groups
  #
  doublezeros <- FALSE
  if (sparse & sm %in% c("RR", "OR") & !(method %in% c("Peto", "GLMM"))) {
    sel.doublezeros <- switch(sm,
                              OR = (event.e == 0   & event.c ==   0) |
                                (event.c == n.c & event.e == n.e),
                              RR = (event.c == 0 & event.e == 0))
    if (any(sel.doublezeros, na.rm = TRUE))
      doublezeros <- TRUE
  }
  #
  # Define continuity correction
  #
  #
  if (avail.incr.both) {
    chknumeric(incr.e, min = 0, NA.ok = FALSE)
    chknumeric(incr.c, min = 0, NA.ok = FALSE)
  }
  else {
    if (addincr) {
      #
      if (is.numeric(incr)) {
        if (is.null(incr.e))
          incr.e <- if (length(incr) == 1) rep(incr, k.all) else incr
        if (is.null(incr.c))
          incr.c <- if (length(incr) == 1) rep(incr, k.all) else incr
      }
      else {
        if (is.tacc) {
          #
          # Treatment arm continuity correction:
          #
          incr.e <- n.e / (n.e + n.c)
          incr.c <- n.c / (n.e + n.c)
        }
      }
    }
    else {
      if (sparse) {
        if (allincr) {
          #
          if (is.numeric(incr)) {
            incr.e <- if (length(incr) == 1) rep(incr, k.all) else incr
            incr.c <- if (length(incr) == 1) rep(incr, k.all) else incr
          }
          else {
            if (is.tacc) {
              #
              # Treatment arm continuity correction:
              #
              incr.e <- n.e / (n.e + n.c)
              incr.c <- n.c / (n.e + n.c)
            }
          }
        }
        else {
          #
          # Bradburn, Deeks, Altman, Stata-procedure "metan":
          # & SAS PROC FREQ (for method = "Inverse")
          #
          if (is.numeric(incr)) {
            incr.e <- incr * sel
            incr.c <- incr * sel
          }
          else {
            if (is.tacc) {
              #
              # Treatment arm continuity correction:
              #
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
    #
    # No continuity correction for Peto method
    #
    if (method == "Peto") {
      incr <- 0
      incr.e <- rep(0, k.all)
      incr.c <- rep(0, k.all)
    }
  }
  #
  incr.e <- incr.e * incl
  incr.c <- incr.c * incl
  #
  if (keepdata) {
    if (missing.subset) {
      data$.incr.e <- incr.e
      data$.incr.c <- incr.c
    }
    else {
      data$.incr.e <- NA
      data$.incr.c <- NA
      #
      data$.incr.e[subset] <- incr.e
      data$.incr.c[subset] <- incr.c
    }
  }
  
  
  ##
  ##
  ## (8) Calculate results for individual studies
  ##
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
    if (method != "Peto") {
      ##
      ## Cooper & Hedges (1994), p. 251-2
      ##
      TE <- log(((n11 + incr.e) * (n22 + incr.c)) /
                ((n12 + incr.e) * (n21 + incr.c)))
      seTE <- sqrt((1 / (n11 + incr.e) + 1 / (n12 + incr.e) +
                    1 / (n21 + incr.c) + 1 / (n22 + incr.c)))
    }
    else {
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
  else if (sm %in% c("RR", "VE")) {
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
  #
  # Set NaN to NA
  #
  TE[is.nan(TE)] <- NA
  
  
  ##
  ##
  ## (9) Additional checks for three-level model
  ##
  ##
  
  three.level <- FALSE
  sel.ni <- !is.infinite(TE) & !is.infinite(seTE)
  ##
  ## Only conduct three-level meta-analysis if variable 'cluster'
  ## contains duplicate values after removing inestimable study
  ## results standard errors
  ##
  if (with.cluster &&
      length(unique(cluster[sel.ni])) != length(cluster[sel.ni]))
    three.level <- TRUE
  ##
  if (three.level) {
    chkmlm(method.tau, missing.method.tau, method.predict,
           method, missing.method)
    ##
    common <- FALSE
    method <- "Inverse"
    is.glmm <- FALSE
    is.lrp <- FALSE
    ##
    if (!(method.tau %in% c("REML", "ML")))
      method.tau <- "REML"
  }
  #
  if (is.lrp) {
    is_installed_package("brglm2", fun, "method", " = \"LRP\"")
    #
    if (!missing.method.tau & method.tau != "DL")
      warn_ignore_input(method.tau, text = "for penalised logistic regression")
    #
    if (by & tau.common)
      stop("Subgroup analysis not defined for penalised logistic regression ",
           "assuming a common tau-squared.",
           call. = FALSE)
  }
  
  
  ##
  ##
  ## (10) Additional checks for GLMM, penalised logistic regression,
  ##      Peto method or SSW
  ##
  ##
  
  if (sm != "OR") {
    if (method == "Peto")
      stop("Peto's method only possible with argument 'sm = \"OR\"'")
    else if (method == "SSW")
      stop("Sample size weighting only available with argument 'sm = \"OR\"'")
    else if (is.glmm)
      stop("Generalised linear mixed models only possible with ",
           "argument 'sm = \"OR\"'.")
    else if (is.lrp)
      stop("Logistic regression with penalised likelihood only possible with ",
           "argument 'sm = \"OR\"'.")
  }
  ##
  if (is.glmm) {
    chkglmm(sm, method.tau, method.random.ci, method.predict,
            adhoc.hakn.ci, adhoc.hakn.pi,
            "OR")
    ##
    if (!is.null(TE.tau)) {
      if (warn)
        warning("Argument 'TE.tau' not considered for GLMM.",
                call. = FALSE)
      TE.tau <- NULL
    }
    ##
    if (!is.null(tau.preset)) {
      if (warn)
        warning("Argument 'tau.preset' not considered for GLMM.",
                call. = FALSE)
      tau.preset <- NULL
    }
    #
    if (model.glmm == "CM.EL")
      is_installed_package("BiasedUrn", fun, "model.glmm", " = \"CM.EL\"")
  }
  #
  if (is.lrp) {
    chklrp(sm, method.tau, method.random.ci, method.predict,
           adhoc.hakn.ci, adhoc.hakn.pi, "OR")
    #
    if (warn) {
      txt.warn <- "for penalised logistic regression"
      #
      warn_ignore_input(TE.tau, !is.null(TE.tau), txt.warn)
      warn_ignore_input(tau.preset, !is.null(tau.preset), txt.warn)
      #
      if (!missing.method.I2 & method.I2 == "tau2")
        warning("Argument 'method.I2' set to \"Q\" for ",
                "penalised logistic regression.",
                call. = FALSE)
    }
    #
    TE.tau <- NULL
    tau.preset <- NULL
    method.I2 <- "Q"
  }
  ##
  ## No need to add anything to cell counts for
  ##  (i)  arcsine difference as summary measure
  ##  (ii) Peto method, GLMM, or penalised logistic regression
  ##
  if (!avail.incr.both) {
    if (sm == "ASD" | method %in% c("Peto", "GLMM", "LRP")) {
      if ((!missing.incr & any(incr != 0)) |
          allincr | addincr |
          (!missing.allstudies & allstudies)
      )
        if (sm == "ASD") {
          if ((sparse | addincr) & warn) {
            warning("Note, no continuity correction considered ",
                    "for arcsine difference (sm = \"ASD\").",
                    call. = FALSE)
          }
        }
      else if (method == "Peto") {
        if ((sparse | addincr) & warn)
          warning("Note, no continuity correction considered ",
                  "for method = \"Peto\".",
                  call. = FALSE)
      }
    }
  }
  
  
  ##
  ##
  ## (11) Do meta-analysis
  ##
  ##
  
  k <- sum(!is.na(event.e[!exclude]) & !is.na(event.c[!exclude]) &
           !is.na(n.e[!exclude]) & !is.na(n.c[!exclude]))
  ##
  for (i in seq_along(method.random.ci))
    if (k == 1 & method.random.ci[i] == "HK")
      method.random.ci[i] <- "classic"
  #
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
      w.common <- C
      TE.common <- log(sum(A, na.rm = TRUE) / sum(C, na.rm = TRUE))
      seTE.common <- sqrt((1 / (2 * sum(A, na.rm = TRUE)^2)  *
                           (sum(A * B, na.rm = TRUE) +
                            exp(TE.common) * (sum(B * C, na.rm = TRUE) +
                                              sum(A * D, na.rm = TRUE)) +
                            exp(TE.common)^2 * sum(C * D, na.rm = TRUE))))
    }
    else if (sm %in% c("RR", "VE")) {
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
      w.common <- S
      TE.common <- log(sum(R, na.rm = TRUE) / sum(S, na.rm = TRUE))
      seTE.common <- sqrt(sum(D, na.rm = TRUE) / (sum(R, na.rm = TRUE) *
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
      w.common <- S
      TE.common <- weighted.mean(TE, w.common, na.rm = TRUE)
      seTE.common <- sqrt(sum(R, na.rm = TRUE) / sum(S, na.rm = TRUE)^2)
    }
    ##
    w.common[is.na(w.common)] <- 0
  }
  else if (method == "Peto") {
    w.common <- 1 / seTE^2
    w.common[exclude] <- 0
    TE.common   <- weighted.mean(TE, w.common, na.rm = TRUE)
    seTE.common <- sqrt(1 / sum(w.common, na.rm = TRUE))
    ##
    w.common[is.na(w.common)] <- 0
  }
  else if (is.glmm) {
    list.bin <- list(ai = event.e[!exclude], n1i = n.e[!exclude],
                     ci = event.c[!exclude], n2i = n.c[!exclude],
                     measure = "OR", model = model.glmm)
    ##
    use.random <-
      sum(!exclude) > 1 &
      !((sum(event.e[!exclude], na.rm = TRUE) == 0 &
         sum(event.c[!exclude], na.rm = TRUE) == 0) |
        (!any(event.e[!exclude] != n.e[!exclude]) |
         !any(event.c[!exclude] != n.c[!exclude])))
    ##
    res.glmm <-
      runGLMM(list.bin,
              method.tau = method.tau,
              method.random.ci = method.random.ci,
              level = level.ma,
              control = list(control),
              use.random = use.random,
              warn = warn)
    ##
    TE.common   <- as.numeric(res.glmm$glmm.common$b)
    seTE.common <- as.numeric(res.glmm$glmm.common$se)
    ##
    w.common <- rep(NA, length(event.e))
  }
  else if (is.lrp) {
    fit.lrp <- runLRP(event.e[!exclude], n1 = n.e[!exclude],
                      event2 = event.c[!exclude], n2 = n.c[!exclude],
                      ...)
    #
    TE.common   <- fit.lrp$TE.common
    seTE.common <- fit.lrp$seTE.common
    #
    w.common <- rep(NA, length(event.e))
  }
  else if (method == "SSW") {
    w.common <- n.e * n.c / (n.e + n.c)
    w.common[exclude] <- 0
    TE.common <- weighted.mean(TE, w.common, na.rm = TRUE)
    seTE.common <- sqrt(sum(w.common^2 * seTE^2, na.rm = TRUE) /
                        sum(w.common, na.rm = TRUE)^2)
    ##
    w.common[is.na(w.common)] <- 0
  }
  ##
  m <- metagen(TE, seTE, studlab,
               exclude = if (missing.exclude) NULL else exclude,
               cluster = cluster, rho = rho,
               #
               weights.common = weights.common,
               weights.random = weights.random,
               #
               sm = sm,
               level = level,
               ##
               common = common,
               random = random,
               overall = overall,
               overall.hetstat = overall.hetstat,
               prediction = prediction,
               ##
               method.tau = if (is.glmm | is.lrp) "DL" else method.tau,
               method.tau.ci = if (is.glmm | is.lrp) "" else method.tau.ci,
               level.hetstat = level.hetstat,
               tau.preset = tau.preset,
               TE.tau = if (Q.Cochrane) TE.common else TE.tau,
               tau.common = FALSE,
               detail.tau = detail.tau,
               #
               method.I2 = method.I2,
               #
               level.ma = level.ma,
               method.common.ci = method.common.ci,
               method.random.ci = method.random.ci,
               adhoc.hakn.ci = adhoc.hakn.ci,
               ##
               level.predict = level.predict,
               method.predict = method.predict,
               adhoc.hakn.pi = adhoc.hakn.pi,
               seed.predict = seed.predict,
               ##
               method.bias = method.bias,
               ##
               backtransf = backtransf,
               ##
               text.common = text.common, text.random = text.random,
               text.predict = text.predict,
               text.w.common = text.w.common, text.w.random = text.w.random,
               ##
               title = title, complab = complab, outclab = outclab,
               #
               label.e = label.e, label.c = label.c,
               label.left = label.left, label.right = label.right,
               col.label.left = col.label.left,
               col.label.right = col.label.right,
               #
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
  #
  # Estimate common tau-squared across subgroups
  #
  if (by & tau.common & !is.glmm)
    hcc <- hetcalc(TE, seTE, method.tau, "",
                   if (Q.Cochrane & method == "MH") TE.common else TE.tau,
                   method.I2, level.hetstat, subgroup, control)
  
  
  ##
  ##
  ## (12) Generate R object
  ##
  ##
  
  res <- list(event.e = event.e, n.e = n.e,
              event.c = event.c, n.c = n.c,
              method = method, method.random = method,
              incr = if (length(unique(incr)) == 1) unique(incr) else incr,
              method.incr = method.incr,
              sparse = sparse,
              allstudies = allstudies,
              doublezeros = doublezeros,
              MH.exact = MH.exact, RR.Cochrane = RR.Cochrane,
              Q.Cochrane = Q.Cochrane,
              Q.CMH = Q.CMH, df.Q.CMH = 1, pval.Q.CMH = pvalQ(Q.CMH, 1),
              print.CMH = print.CMH,
              incr.e = incr.e, incr.c = incr.c,
              k.MH = if (method == "MH") sum(w.common > 0) else NA)
  ##
  ## Add meta-analysis results
  ## (after removing unneeded list elements)
  ##
  m$method <- NULL
  m$method.random <- NULL
  m$n.e <- NULL
  m$n.c <- NULL
  m$pscale <- NULL
  m$irscale <- NULL
  m$irunit <- NULL
  m$method.ci <- NULL
  m$method.mean <- NULL
  m$approx.TE <- NULL
  m$approx.seTE <- NULL
  ##
  res <- c(res, m)
  ##
  ## Add data
  ##
  res$TE.tau <- TE.tau
  ##
  res$pscale <- pscale
  #
  if (is.pairwise | data.pairwise) {
    res$pairwise <- TRUE
    res$k.study <- length(unique(res$studlab[!is.na(res$TE)]))
  }
  #
  res$call <- match.call()
  ##
  if (method %in% c("MH", "Peto", "GLMM", "LRP", "SSW")) {
    res <- ci2meta(res, ci.c = ci(TE.common, seTE.common, level = level.ma))
    res$w.common <- w.common
  }
  #
  if (is.glmm) {
    res$method.tau <- method.tau
    res <- addGLMM(res, res.glmm, method.I2)
    res$model.glmm <- model.glmm
    ##
    if (by) {
      n.subgroups <- length(unique(subgroup[!exclude]))
      if (n.subgroups > 1)
        subgroup.glmm <-
          factor(subgroup[!exclude], bylevs(subgroup[!exclude]))
      ##
      hcc <-
        hccGLMM(
          res,
          runGLMM(list.bin,
                  method.tau = method.tau,
                  method.random.ci = method.random.ci,
                  level = level.hetstat,
                  data =
                    if (n.subgroups > 1)
                      list(data.frame(subgroup.glmm))
                    else NULL,
                  mods =
                    if (n.subgroups > 1)
                      as.call(~ subgroup.glmm)
                    else
                      NULL,
                  control = list(control),
                  use.random = use.random,
                  warn = warn)$glmm.random[[1]],
          method.I2
        )
    }
  }
  else if (is.lrp) {
    res <- ci2meta(res,
                   ci.r = ci(fit.lrp$TE.random, fit.lrp$seTE.random,
                             level = level.ma,
                             df = ifelse(method.random.ci == "HK",
                                         m$k - 1, Inf)))
    res$w.random <- w.common
    res$phi <- fit.lrp$phi
  }
  else if (method == "SSW") {
    res <- ci2meta(res,
                   ci.r = ci(TE.random, seTE.random, level = level.ma,
                             df = ifelse(method.random.ci == "HK",
                                         m$k - 1, Inf)))
    res$w.random <- w.random
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
    res$subgroup <- subgroup
    res$subgroup.name <- subgroup.name
    res$print.subgroup.name <- print.subgroup.name
    res$sep.subgroup <- sep.subgroup
    res$test.subgroup <- test.subgroup
    res$prediction.subgroup <- prediction.subgroup
    res$tau.common <- tau.common
    ##
    if (!tau.common) {
      res <- c(res, subgroup(res, seed = seed.predict.subgroup))
      if (res$three.level)
        res <- setNA3(res)
    }
    else if (!is.null(tau.preset))
      res <-
        c(res, subgroup(res, tau.preset, seed = seed.predict.subgroup))
    else {
      if (is.glmm)
        res <- c(res,
                 subgroup(res, NULL,
                          factor(res$subgroup, bylevs(res$subgroup)), ...))
      else if (res$three.level)
        res <- c(res,
                 subgroup(res, NULL,
                          factor(res$subgroup, bylevs(res$subgroup))))
      else
        res <-
          c(res, subgroup(res, hcc$tau.resid, seed = seed.predict.subgroup))
    }
    ##
    if (tau.common && is.null(tau.preset))
      res <- addHet(res, hcc, !(is.glmm | is.lrp))
    ##
    res$n.w <- NULL
    res$event.w <- NULL
    ##
    res$n.harmonic.mean.w <- NULL
    ##
    res$time.e.w <- NULL
    res$time.c.w <- NULL
    res$t.harmonic.mean.w <- NULL
    ##
    res <- setNAwithin(res, res$three.level | is.glmm)
  }
  #
  # Mantel-Haenszel method is common effect method
  #
  if (res$method.random == "MH")
    res$method.random <- "Inverse"
  #
  # Do not return tau^2 and tau for penalised logistic regression
  #
  if (is.lrp) {
    res$method.random <- "LRP"
    #
    res$tau2 <- res$lower.tau2 <- res$upper.tau2 <- NA
    res$tau <- res$lower.tau <- res$upper.tau <- NA
    res$method.tau <- ""
    #
    res <- calcPI(res)
    #
    res$version.brglm2 <- packageDescription("brglm2")$Version
  }
  ##
  ## Backward compatibility
  ##
  res <- backward(res)
  ##
  class(res) <- c(fun, "meta")
  
  
  res
}
