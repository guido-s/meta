#' Mental Health Treatment
#' 
#' @description
#' Meta-analysis on the Effect of Mental Health Treatment on Medical
#' Utilisation.
#' 
#' Data example in Fleiss (1993) for meta-analysis with continuous outcomes.
#' 
#' @name Fleiss93cont
#' 
#' @docType data
#' 
#' @format A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study label \cr
#' \bold{\emph{year}}\tab year of publication \cr
#' \bold{\emph{n.e}}\tab number of observations in psychotherapy group
#'   \cr
#' \bold{\emph{mean.e}}\tab estimated mean in psychotherapy group \cr
#' \bold{\emph{sd.e}}\tab standard deviation in psychotherapy group
#'   \cr
#' \bold{\emph{n.c}}\tab number of observations in control group \cr
#' \bold{\emph{mean.c}}\tab estimated mean in control group \cr
#' \bold{\emph{sd.c}}\tab standard deviation in control group
#' }
#' 
#' @seealso \code{\link{Fleiss93}}
#' 
#' @source
#' Fleiss JL (1993):
#' The statistical basis of meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{2}, 121--45
#'
#' @keywords datasets
#'
#' @examples
#' data(Fleiss93cont)
#' metacont(n.e, mean.e, sd.e,
#'          n.c, mean.c, sd.c,
#'          data = Fleiss93cont,
#'          studlab = paste(study, year),
#'          comb.random = FALSE)


NULL
