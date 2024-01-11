#' Mental Health Treatment
#' 
#' @description
#' Meta-analysis on the Effect of Mental Health Treatment on Medical
#' Utilisation.
#' 
#' Data example in Fleiss (1993) for meta-analysis with continuous outcomes.
#' 
#' @name Fleiss1993cont
#' @aliases Fleiss93cont
#' 
#' @docType data
#' 
#' @format A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study label \cr
#' \bold{\emph{year}}\tab year of publication \cr
#' \bold{\emph{n.psyc}}\tab number of observations in psychotherapy
#'   group \cr
#' \bold{\emph{mean.psyc}}\tab estimated mean in psychotherapy group
#'   \cr
#' \bold{\emph{sd.psyc}}\tab standard deviation in psychotherapy group
#'   \cr
#' \bold{\emph{n.cont}}\tab number of observations in control group
#'   \cr
#' \bold{\emph{mean.cont}}\tab estimated mean in control group \cr
#' \bold{\emph{sd.cont}}\tab standard deviation in control group
#' }
#' 
#' @seealso \code{\link{Fleiss1993bin}}
#' 
#' @source
#' Fleiss JL (1993):
#' The statistical basis of meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{2}, 121--45
#'
#' @keywords datasets
#' 
#' @seealso \code{\link{metacont}}
#'
#' @examples
#' data(Fleiss1993cont)
#' # Note, the following command uses the bias-corrected version of
#' # Hedges' g. Accordingly, results differ from Fleiss (1993), section 3,
#' # using the uncorrected version of Hedges' g.
#' metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, studlab = paste(study, year),
#'   random = FALSE, sm = "SMD")


NULL
