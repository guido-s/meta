#' Description of summary measures available in R package \bold{meta}
#' 
#' @description
#' Description of summary measures available in R package \bold{meta}
#' 
#' @details
#' The following summary measures (argument \code{sm}) are recognised
#' in R package \bold{meta}.
#'
#' \subsection{Meta-analysis of binary outcome data (\code{\link{metabin})}}{
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Summary measure} \cr
#' \code{sm = "OR"} \tab Odds ratio (Fleiss, 1993)\cr
#' \code{sm = "RR"} \tab Risk ratio (Fleiss, 1993) \cr
#' \code{sm = "RD"} \tab Risk difference (Fleiss, 1993) \cr
#' \code{sm = "ASD"} \tab Arcsine difference (Rücker et al., 2009) \cr
#' \code{sm = "DOR"} \tab Diagnostic odds ratio (Moses et al., 1993)
#'   \cr
#' \code{sm = "VE"} \tab Vaccine efficacy or vaccine effectiveness
#' }
#'
#' Note, mathematically, odds ratios and diagnostic odds ratios are
#' identical, however, the labels in printouts and figures
#' differ. Furthermore, log risk ratio (logRR) and log vaccine ratio
#' (logVR) are mathematical identical, however, back-transformed
#' results differ as vaccine efficacy or effectiveness is defined as
#' \code{VE = 100 * (1 - RR)}.
#'
#' A continuity correction is used for some summary measures in the
#' case of a zero cell count (see \code{\link{metabin}}). 
#'
#' List elements \code{TE}, \code{TE.common}, \code{TE.random}, etc.,
#' contain transformed values, e.g., log odds ratios, log risk ratios
#' or log vaccine ratios.  In printouts and plots transformed values
#' are back transformed if argument \code{backtransf = TRUE}
#' (default), with exception of the arcsine difference where no
#' back-transformation exists. Auxiliary function
#' \code{\link{logVR2VE}} is used to back-transform log vaccine ratios
#' to vaccine efficacy or effectiveness while \code{\link[base]{exp}}
#' is used to back-transform log odds or risk ratios.
#' }
#'
#' \subsection{Meta-analysis of continuous outcome data (\code{\link{metacont})}}{
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Summary measure} \cr
#' \code{sm = "MD"} \tab Mean difference \cr
#' \code{sm = "SMD"} \tab Standardised mean difference \cr
#' \code{sm = "ROM"} \tab Ratio of means
#' }
#'
#' Three variants to calculate the standardised mean difference are
#' available (see \code{\link{metacont}}).
#' 
#' For the ratio of means, list elements \code{TE}, \code{TE.common},
#' \code{TE.random}, etc., contain the log transformed ratio of
#' means. In printouts and plots these values are back transformed
#' using \code{\link[base]{exp}} if argument \code{backtransf = TRUE}
#' (default).
#' }
#'
#' \subsection{Meta-analysis of correlations (\code{\link{metacor})}}{
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Summary measure} \cr
#' \code{sm = "ZCOR"} \tab Fisher's z transformed correlation \cr
#' \code{sm = "COR"} \tab Untransformed correlations
#' }
#'
#' For Fisher's z transformed correlations, list elements \code{TE},
#' \code{TE.common}, \code{TE.random}, etc., contain the transformed
#' correlations. In printouts and plots these values are back
#' transformed using auxiliary function \code{\link{z2cor}} if
#' argument \code{backtransf = TRUE} (default).
#' }
#'
#' \subsection{Meta-analysis of incidence rates (\code{\link{metainc})}}{
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Summary measure} \cr
#' \code{sm = "IRR"} \tab Incidence rate ratio \cr
#' \code{sm = "IRD"} \tab Incidence rate difference \cr
#' \code{sm = "IRSD"} \tab Square root transformed incidence rate
#'   difference \cr
#' \code{sm = "VE"} \tab Vaccine efficacy or vaccine effectiveness
#' }
#'
#' Note, log incidence rate ratio (logIRR) and log vaccine ratio
#' (logVR) are mathematical identical, however, back-transformed
#' results differ as vaccine efficacy or effectiveness is defined as
#' \code{VE = 100 * (1 - IRR)}.
#'
#' List elements \code{TE}, \code{TE.common}, \code{TE.random}, etc.,
#' contain the transformed incidence rates. In printouts and plots
#' these values are back transformed if argument \code{backtransf =
#' TRUE} (default). For back-transformation, \code{\link[base]{exp}}
#' is used for the incidence rate ratio, power of 2 is used for square
#' root transformed rates and \code{\link{logVR2VE}} is used for
#' vaccine efficacy / effectiveness.
#' }
#'
#' \subsection{Meta-analysis of single means (\code{\link{metamean})}}{
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Summary measure} \cr
#' \code{sm = "MRAW"} \tab Raw, i.e. untransformed, means \cr
#' \code{sm = "MLN"} \tab Log transformed means
#' }
#' 
#' Calculations are conducted on the log scale if \code{sm =
#' "MLN"}. Accordingly, list elements \code{TE}, \code{TE.common}, and
#' \code{TE.random} contain the logarithm of means. In printouts and
#' plots these values are back transformed using
#' \code{\link[base]{exp}} if argument \code{backtransf = TRUE}.
#' }
#'
#' \subsection{Meta-analysis of single proportions (\code{\link{metaprop})}}{
#' 
#' The following transformations of proportions are
#' implemented to calculate an overall proportion:
#' 
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Summary measure} \cr
#' \code{sm = "PLOGIT"} \tab Logit transformation \cr
#' \code{sm = "PAS"} \tab Arcsine transformation \cr
#' \code{sm = "PFT"} \tab Freeman-Tukey Double arcsine transformation
#'   \cr
#' \code{sm = "PLN"} \tab Log transformation \cr
#' \code{sm = "PRAW"} \tab No transformation
#' }
#'
#' List elements \code{TE}, \code{TE.common}, \code{TE.random}, etc.,
#' contain the transformed proportions. In printouts and plots these
#' values are back transformed if argument \code{backtransf = TRUE}
#' (default). For back-transformation, \code{\link{logit2p}} is used
#' for logit transformed proportions, \code{\link{asin2p}} is used for
#' (Freeman-Tukey) arcsine transformed proportions and
#' \code{\link[base]{exp}} is used for log transformed proportions.
#' }
#'
#' \subsection{Meta-analysis of single rates (\code{\link{metarate})}}{
#' 
#' The following transformations of incidence rates are implemented to
#' calculate an overall rate:
#' 
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Summary measure} \cr
#' \code{sm = "IRLN"} \tab Log transformation \cr
#' \code{sm = "IRS"} \tab Square root transformation \cr
#' \code{sm = "IRFT"} \tab Freeman-Tukey Double arcsine transformation
#'   \cr
#' \code{sm = "IR"} \tab No transformation
#' }
#'
#' List elements \code{TE}, \code{TE.common}, \code{TE.random}, etc.,
#' contain the transformed incidence rates. In printouts and plots
#' these values are back transformed if argument \code{backtransf =
#' TRUE} (default). For back-transformation, \code{\link[base]{exp}}
#' is used for log transformed rates, power of 2 is used for square
#' root transformed rates and \code{\link{asin2ir}} is used for
#' Freeman-Tukey arcsine transformed rates.
#' }
#'
#' \subsection{Generic inverse variance method (\code{\link{metagen})}}{
#' 
#' The following summary measures are recognised in addition to the
#' above mentioned summary measures:
#' 
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Summary measure} \cr
#' \code{sm = "HR"} \tab Hazard ratio \cr
#' \code{sm = "VE"} \tab Vaccine efficacy or vaccine effectiveness
#' }
#'
#' List elements \code{TE}, \code{TE.common}, \code{TE.random}, etc.,
#' contain transformed values, i.e., log hazard ratios and log vaccine
#' ratios. In printouts and plots these values are back transformed if
#' argument \code{backtransf = TRUE} (default); see also
#' \link{meta-transf}.
#' }
#'
#' @name meta-sm
#' 
#' @aliases meta-sm
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{meta-package}}, \code{\link{meta-transf}},
#'   \code{\link{meta-object}}, \code{\link{print.meta}},
#'   \code{\link{summary.meta}}, \code{\link{forest.meta}}
#' 
#' @references
#' Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2010):
#' A basic introduction to fixed-effect and random-effects models for
#' meta-analysis.
#' \emph{Research Synthesis Methods},
#' \bold{1}, 97--111
#' 
#' Fleiss JL (1993):
#' The statistical basis of meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{2}, 121--45
#'
#' Moses LE, Shapiro D, Littenberg B (1993):
#' Combining Independent Studies of a Diagnostic Test into a Summary
#' Roc Curve: Data-Analytic Approaches and Some Additional
#' Considerations.
#' \emph{Statistics in Medicine},
#' \bold{12}, 1293--1316
#' 
#' Rücker G, Schwarzer G, Carpenter J, Olkin I (2009):
#' Why add anything to nothing? The arcsine difference as a measure of
#' treatment effect in meta-analysis with zero cells.
#' \emph{Statistics in Medicine},
#' \bold{28}, 721--38
#' 
#' Stijnen T, Hamza TH, Ozdemir P (2010):
#' Random effects meta-analysis of event outcome in the framework of
#' the generalized linear mixed model with applications in sparse
#' data.
#' \emph{Statistics in Medicine},
#' \bold{29}, 3046--67


NULL
