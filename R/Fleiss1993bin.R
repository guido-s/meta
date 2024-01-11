#' Aspirin after Myocardial Infarction
#' 
#' @description
#' Meta-analysis on aspirin in preventing death after myocardial
#' infarction.
#' 
#' Data example in Fleiss (1993) for meta-analysis with binary
#' outcomes.
#' 
#' @name Fleiss1993bin
#' @aliases Fleiss93
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study label \cr
#' \bold{\emph{year}}\tab year of publication \cr
#' \bold{\emph{d.asp}}\tab number of deaths in aspirin group \cr
#' \bold{\emph{n.asp}}\tab number of observations in aspirin group \cr
#' \bold{\emph{d.plac}}\tab number of deaths in placebo group \cr
#' \bold{\emph{n.plac}}\tab number of observations in placebo group
#' }
#' 
#' @source
#' Fleiss JL (1993):
#' The statistical basis of meta-analysis.
#' \emph{Statistical Methods in Medical Research},
#' \bold{2}, 121--45
#' 
#' @keywords datasets
#' 
#' @seealso \code{\link{metabin}}
#' 
#' @examples
#' data(Fleiss1993bin)
#' metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin,
#'   studlab = paste(study, year), sm = "OR", random = FALSE)


NULL
