#' Meta-analysis on Prevention of First Bleeding in Cirrhosis
#' 
#' @description
#' Meta-analysis on Prevention of First Bleeding in Cirrhosis
#' comparing beta-blocker or sclerotherapy with placebo.
#' 
#' @name Pagliaro1992
#' @aliases Pagliaro1992
#' 
#' @docType data
#' 
#' @format A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{id}}\tab study id \cr
#' \bold{\emph{treat.exp}}\tab treatment in experimental group \cr
#' \bold{\emph{logOR}}\tab log odds ratio
#'   \cr
#' \bold{\emph{selogOR}}\tab standard error of log odds ratio
#'   \cr
#' \bold{\emph{bleed.exp}}\tab number of bleedings in experimental group \cr
#' \bold{\emph{n.cont}}\tab number of observations in experimental group \cr
#' \bold{\emph{bleed.plac}}\tab number of bleedings in placebo group \cr
#' \bold{\emph{n.plac}}\tab number of observations in placebo group
#' }
#' 
#' @source
#' Pagliaro L, D’Amico G et al. (1992):
#' Prevention of first bleeding in cirrhosis.
#' \emph{Annals in Internal Medicine},
#' \bold{117}, 59--70 
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Pagliaro1992)
#' sclero <- subset(Pagliaro1992, treat.exp == "Sclerotherapy")
#' 
#' m <- metagen(logOR, selogOR, data = sclero, sm = "OR")
#' m
#'
#' # Thompson & Sharp (1999), Table IV, method (2)
#' metabias(m, method = "Egger")
#'
#' # Thompson & Sharp (1999), Table IV, method (3a)
#' metabias(m, method = "Thompson")
#'
#' # Thompson & Sharp (1999), Table IV, method (3b)
#' update(m, method.tau = "ML")
#' metabias(update(m, method.tau = "ML"), method = "Thompson")
#' 


NULL
