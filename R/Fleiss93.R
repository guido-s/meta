#' Aspirin after Myocardial Infarction
#' 
#' @description
#' Meta-analysis on aspirin in preventing death after myocardial
#' infarction.
#' 
#' Data example in Fleiss (1993) for meta-analysis with binary
#' outcomes.
#' 
#' @name Fleiss93
#' 
#' @docType data
#' 
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study label \cr
#' \bold{\emph{year}}\tab year of publication \cr
#' \bold{\emph{event.e}}\tab number of deaths in aspirin group \cr
#' \bold{\emph{n.e}}\tab number of observations in aspirin group \cr
#' \bold{\emph{event.c}}\tab number of deaths in placebo group \cr
#' \bold{\emph{n.c}}\tab number of observations in placebo group
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
#' @examples
#' data(Fleiss93)
#' metabin(event.e, n.e, event.c, n.c,
#'         data = Fleiss93,
#'         studlab = paste(study, year),
#'         sm = "OR", comb.random = FALSE)


NULL
