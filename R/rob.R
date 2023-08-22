#' Risk of bias assessment
#' 
#' @description
#' Create table with risk of bias assessment or add table to existing
#' meta-analysis
#' 
#' @param item1 Risk of bias item 1.
#' @param item2 Risk of bias item 2.
#' @param item3 Risk of bias item 3.
#' @param item4 Risk of bias item 4.
#' @param item5 Risk of bias item 5.
#' @param item6 Risk of bias item 6.
#' @param item7 Risk of bias item 7.
#' @param item8 Risk of bias item 8.
#' @param item9 Risk of bias item 9.
#' @param item10 Risk of bias item 10.
#' @param studlab Study labels.
#' @param data A data frame or an object of class \code{meta}.
#' @param method Risk of bias (RoB) tool.
#' @param domains A character vector with names of RoB domains.
#' @param categories Possible RoB categories.
#' @param col Colours for RoB categories.
#' @param symbols Corresponding symbols for RoB categories.
#' @param legend A logical specifying whether legend with RoB domains
#'   should be printed.
#' @param overwrite A logical indicating whether an existing risk of
#'   bias table in a meta-analysis object should be overwritten.
#' @param x An object of class \code{rob}.
#' @param \dots Additional printing arguments
#' 
#' @details
#' This function can be used to define a risk of bias (RoB) assessment
#' for a meta-analysis which can be shown in a forest plot. The
#' resulting risk of bias table contains study labels and variables
#' for the RoB domains; variable names for RoB domains are equal to A,
#' B, etc. The RoB table is directly returned if argument \code{data}
#' is a data frame. The RoB table is added as a new list element 'rob'
#' to a meta-analysis object if argument \code{data} is a
#' meta-analysis.
#'
#' The user must either specify the categories and (optionally)
#' domains of the RoB tool (using the eponymous arguments) or one of
#' the following RoB tools.
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Risk of bias tool} \cr
#' \code{method = "RoB1"} \tab RoB 1 tool for randomized studies \cr
#' \code{method = "RoB2"} \tab RoB 2 tool for randomized studies \cr
#' \code{method = "RoB2-cluster"} \tab RoB 2 tool for
#'   cluster-randomized trials \cr
#' \code{method = "RoB2-crossover"} \tab RoB 2 tool for crossover
#'   trials \cr
#' \code{method = "ROBINS-I"} \tab Risk Of Bias In Non-randomized
#'   Studies - of Interventions \cr
#' \code{method = "ROBINS-E"} \tab Risk Of Bias In Non-randomized
#'   Studies - of Exposures
#' }
#' These RoB tools are described on the website
#' \url{https://www.riskofbias.info/}.
#'
#' \subsection{Risk of bias domains}{
#'
#' By default, i.e., if argument \code{domains} is not provided by the
#' user, the following names are used for RoB domains.
#'
#' \itemize{
#' \item RoB 1 tool for randomized studies (RoB1):
#'
#' \enumerate{
#'  \item Random sequence generation (selection bias)
#'  \item Allocation concealment (selection bias)
#'  \item Blinding of participants and personnel (performance bias)
#'  \item Incomplete outcome data (attrition bias)
#'  \item Selective reporting (reporting bias)
#'  \item Other bias
#' }
#'
#' \item RoB 2 tool for randomized studies (RoB2):
#' \enumerate{
#'  \item Bias arising from the randomization process
#'  \item Bias due to deviations from intended intervention"
#'  \item Bias due to missing outcome data
#'  \item Bias in measurement of the outcome
#'  \item Bias in selection of the reported result
#' }
#'
#' \item RoB 2 tool for cluster-randomized trials (RoB2-cluster):
#' \enumerate{
#'  \item Bias arising from the randomization process
#'  \item Bias arising from the identification or recruitment of
#'    participants into clusters
#'  \item Bias due to deviations from intended intervention
#'  \item Bias due to missing outcome data
#'  \item Bias in measurement of the outcome
#'  \item Bias in selection of the reported result
#' }
#'
#' \item RoB 2 tool for crossover trials (RoB2-crossover)
#' \enumerate{
#'  \item Bias arising from the randomization process
#'  \item Bias arising from period and carryover effects
#'  \item Bias due to deviations from intended intervention
#'  \item Bias due to missing outcome data
#'  \item Bias in measurement of the outcome
#'  \item Bias in selection of the reported result
#' }
#' 
#' \item Risk Of Bias In Non-randomized Studies - of Intervention (ROBINS-I):
#' \enumerate{
#'  \item Risk of bias due to confounding
#'  \item Risk of bias in selection of participants into the study
#'  \item Risk of bias in classification of interventions
#'  \item Risk of bias due to deviations from intented interventions
#'  \item Risk of bias due to missing outcome data
#'  \item Risk of bias in measurement of the outcome
#'  \item Risk of bias in the selection of the reported results
#' }
#' 
#' \item Risk Of Bias In Non-randomized Studies - of Exposures (ROBINS-E):
#' \enumerate{
#'  \item Risk of bias due to confounding
#'  \item Risk of bias arising from measurement of the exposure into
#'    the study (or into the analysis)
#'  \item Risk of bias due to post-exposure interventions
#'  \item Risk of bias due to deviations from intented interventions
#'  \item Risk of bias due to missing outcome data
#'  \item Risk of bias in measurement of the outcome
#'  \item Risk of bias in the selection of the reported results
#' }
#' 
#' \item User-defined RoB assessment:
#' \enumerate{
#'  \item Item 1
#'  \item Item 2
#'  \item \dots
#' }
#' }
#'
#' It is possible to define additional bias domains for the available
#' RoB tools. In this case, only the names for new RoB domains have to
#' be provided in argument \code{domains}. If argument \code{domains}
#' is not used to specify new domains, the names "Additional item 1"
#' etc. will be used. It is also possible to modify the pre-defined
#' domain names using argument \code{domains}.
#'
#' 
#' The maximum number of bias domains / items is ten (see arguments
#' \code{item1}, ..., \code{item10}).
#' }
#'
#' \subsection{Risk of bias categories, colours and symbols}{
#' 
#' By default, the following settings are used.
#'
#' RoB 1 tool:
#'
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Values} \cr
#' \code{categories} \tab "Low risk of bias", "Unclear risk of bias",
#'   "High risk of bias" \cr
#' \code{col} \tab "green", "yellow", "red" \cr
#' \code{symbols} \tab "+", "?", "-"
#' }
#'
#' RoB 2 tools:
#'
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Values} \cr
#' \code{categories} \tab "Low risk of bias", "Some concerns",
#'   "High risk of bias" \cr
#' \code{col} \tab "green", "yellow", "red" \cr
#' \code{symbols} \tab "+", "?", "-"
#' }
#'
#' ROBINS tools:
#'
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Values} \cr
#' \code{categories} \tab "Low risk", "Some concerns", "High risk",
#'   "Very high risk", "NI" \cr
#' \code{col} \tab "green", "yellow", "red", "darkred", "darkgrey" \cr
#' \code{symbols} \tab none
#' }
#'
#' User-defined RoB tools:
#'
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Values} \cr
#' \code{categories} \tab Must be specified by the user \cr
#' \code{col} \tab 1, 2, ... \cr
#' \code{symbols} \tab none
#' }
#'
#' If colours (\code{col}) and symbols (\code{symbols}) are provided,
#' they must be of same length as the number of categories. 
#' }
#' 
#' @return
#' A data frame with study labels and risk of bias items and
#' additional class "rob".
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.meta}}, \code{\link{metagen}}
#' 
#' @examples
#' data(Fleiss1993cont)
#' # Do meta-analysis (common effect and random effects model)
#' #
#' meta1 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, studlab = paste(study, year), sm = "SMD")
#' 
#' # A fictive risk of bias assessment with three domains
#' #
#' meta2 <- rob(n.psyc > 20, mean.psyc > 5, sd.psyc > 2,
#'   data = meta1, categories = c(FALSE, TRUE), col = c("green", "red"))
#' forest(meta2, rob.lab = "RoB")
#' 
#' \dontrun{
#' # Risk of bias 1 tool (with fictive RoB assessment)
#' set.seed(1909)
#' Fleiss1993cont$D1 <-
#'   sample(c("low", "unclear", "high"), meta1$k, replace = TRUE)
#' Fleiss1993cont$D2 <-
#'   sample(c("low", "unclear", "high"), meta1$k, replace = TRUE)
#' Fleiss1993cont$D3 <-
#'   sample(c("low", "unclear", "high"), meta1$k, replace = TRUE)
#' Fleiss1993cont$D4 <-
#'   sample(c("low", "unclear", "high"), meta1$k, replace = TRUE)
#' Fleiss1993cont$D5 <-
#'   sample(c("low", "unclear", "high"), meta1$k, replace = TRUE)
#' Fleiss1993cont$D6 <-
#'   sample(c("low", "unclear", "high"), meta1$k, replace = TRUE)
#' meta3 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
#'   data = Fleiss1993cont, studlab = paste(study, year), sm = "SMD")
#' meta4 <- rob(D1, D2, D3, D4, D5, D6, data = meta3, method = "rob1")
#' forest(meta4)
#' }
#' 
#' @export rob


rob <- function(item1,
                item2 = NULL,
                item3 = NULL,
                item4 = NULL,
                item5 = NULL,
                item6 = NULL,
                item7 = NULL,
                item8 = NULL,
                item9 = NULL,
                item10 = NULL,
                studlab = NULL,
                data = NULL,
                ##
                method = gs("method.rob"),
                domains = NULL,
                categories = NULL,
                col = NULL,
                symbols = NULL,
                legend = TRUE,
                ##
                overwrite = FALSE) {
  
  ##
  ##
  ## (1) Read data
  ##
  ##
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  meta.object <- NULL
  ##
  is.meta <- inherits(data, "meta")
  ##
  if (nulldata)
    data <- sfsp
  else if (is.meta) {
    meta.object <- data
    data <- data$data
    if (!is.null(meta.object$subset))
      data <- data[meta.object$subset, ]
    if (is.null(studlab))
      studlab <- data$.studlab
  }
  ##
  ## Catch 'item1', etc.
  ##
  item1 <- catch("item1", mc, data, sfsp)
  chknull(item1)
  k.All <- length(item1)
  ##
  item2 <- catch("item2", mc, data, sfsp)
  item3 <- catch("item3", mc, data, sfsp)
  item4 <- catch("item4", mc, data, sfsp)
  item5 <- catch("item5", mc, data, sfsp)
  item6 <- catch("item6", mc, data, sfsp)
  item7 <- catch("item7", mc, data, sfsp)
  item8 <- catch("item8", mc, data, sfsp)
  item9 <- catch("item9", mc, data, sfsp)
  item10 <- catch("item10", mc, data, sfsp)
  ##
  studlab <- catch("studlab", mc, data, sfsp)
  
  
  ##
  ##
  ## (2) Check data
  ##
  ##
  
  fun <- "rob"
  ##
  if (!is.null(item2))
    chklength(item2, k.All, fun)
  ##
  if (!is.null(item3)) {
    if (is.null(item2))
      stop("Argument 'item3' provided, but argument 'item2' is NULL.",
           call. = FALSE)
    chklength(item3, k.All, fun)
  }
  ##
  if (!is.null(item4)) {
    if (is.null(item3))
      stop("Argument 'item4' provided, but argument 'item3' is NULL.",
           call. = FALSE)
    chklength(item4, k.All, fun)
  }
  ##
  if (!is.null(item5)) {
    if (is.null(item4))
      stop("Argument 'item5' provided, but argument 'item4' is NULL.",
           call. = FALSE)
    chklength(item5, k.All, fun)
  }
  ##
  if (!is.null(item6)) {
    if (is.null(item5))
      stop("Argument 'item6' provided, but argument 'item5' is NULL.",
           call. = FALSE)
    chklength(item6, k.All, fun)
  }
  ##
  if (!is.null(item7)) {
    if (is.null(item6))
      stop("Argument 'item7' provided, but argument 'item6' is NULL.",
           call. = FALSE)
    chklength(item7, k.All, fun)
  }
  ##
  if (!is.null(item8)) {
    if (is.null(item7))
      stop("Argument 'item8' provided, but argument 'item7' is NULL.",
           call. = FALSE)
    chklength(item8, k.All, fun)
  }
  ##
  if (!is.null(item9)) {
    if (is.null(item8))
      stop("Argument 'item9' provided, but argument 'item8' is NULL.",
           call. = FALSE)
    chklength(item9, k.All, fun)
  }
  ##
  if (!is.null(item10)) {
    if (is.null(item9))
      stop("Argument 'item10' provided, but argument 'item9' is NULL.",
           call. = FALSE)
    chklength(item10, k.All, fun)
  }
  ##
  if (is.null(studlab)) {
    if (is.meta)
      studlab <- data$.studlab
    else
      studlab <- seq_len(k.All)
  }
  else
    chklength(studlab, k.All, fun)
  ##
  chklogical(legend)
  chklogical(overwrite)
  ##
  if (!is.null(method))
    method <- setchar(method, gs("meth4rob"))
  else
    method <- "user-defined"
  ##
  is.RoB <- substring(method, 1, 3) == "RoB"
  is.ROBINS <- substring(method, 1, 6) == "ROBINS"
  ##
  if (method == "RoB1")
    dm <-
      c("(A) Random sequence generation (selection bias)",
        "(B) Allocation concealment (selection bias)",
        "(C) Blinding of participants and personnel (performance bias)",
        "(D) Incomplete outcome data (attrition bias)",
        "(E) Selective reporting (reporting bias)",
        "(F) Other bias")
  else if (method == "RoB2")
    dm <-
      c("(A) Bias arising from the randomization process",
        "(B) Bias due to deviations from intended intervention",
        "(C) Bias due to missing outcome data",
        "(D) Bias in measurement of the outcome",
        "(E) Bias in selection of the reported result")
  else if (method == "RoB2-cluster")
    dm <-
      c("(A) Bias arising from the randomization process",
        paste("(B) Bias arising from the identification or",
              "recruitment of participants into clusters"),
        "(C) Bias due to deviations from intended intervention",
        "(D) Bias due to missing outcome data",
        "(E) Bias in measurement of the outcome",
        "(F) Bias in selection of the reported result")
  else if (method == "RoB2-crossover")
    dm <-
      c("(A) Bias arising from the randomization process",
        "(B) Bias arising from period and carryover effects",
        "(C) Bias due to deviations from intended intervention",
        "(D) Bias due to missing outcome data",
        "(E) Bias in measurement of the outcome",
        "(F) Bias in selection of the reported result")
  else if (method == "ROBINS-I")
    dm <-
      c("(A) Risk of bias due to confounding",
        "(B) Risk of bias in selection of participants into the study",
        "(C) Risk of bias in classification of interventions",
        "(D) Risk of bias due to deviations from intented interventions",
        "(E) Risk of bias due to missing outcome data",
        "(F) Risk of bias in measurement of the outcome",
        "(G) Risk of bias in the selection of the reported results")
  else if (method == "ROBINS-E")
    dm <-
      c("(A) Risk of bias due to confounding",
        "(B) Risk of bias arising from measurement of the exposure",
        paste("(C) Risk of bias in selection of participants",
              "into the study (or into the analysis)"),
        "(D) Risk of bias due to post-exposure interventions",
        "(E) Risk of bias due to missing data",
        "(F) Risk of bias in measurement of the outcome",
        "(G) Risk of bias in selection of the reported results")
  else
    dm <- NULL
  
  
  ##
  ##
  ## (3) Internal functions
  ##
  ##
  
  setcat <- function(x, labels) {
    if (is.null(labels))
      return(x)
    ##
    x <- setchar(x, labels)
    factor(x, levels = labels)
  }
  ##
  setdom <- function(x, method, dm, n.cols) {
    
    n.domains <- length(x)
    ##
    if (method == "user-defined") {
      if (n.domains != n.cols)
        stop("Number of domain names does not match number of domains.",
             call. = FALSE)
      else
        return(x)
    }
    ##
    if (method %in% c("RoB1", "RoB2-cluster", "RoB2-crossover")) {
      if (n.domains != n.cols) {
        if (n.domains == n.cols - 6) {
          x <- paste0("(", LETTERS[6 + seq_len(n.cols - 6)], ") ", x)
          return(c(dm, x))
        }
        else if (n.domains < n.cols) {
          x <- paste0("(", LETTERS[6 + seq_len(n.cols - 6)], ") ",
                      "Additional item ", seq_len(n.cols - 6))
          return(c(dm, x))
        }
        else
          stop("Wrong number of domains provided for '", method,
               "' (must be ", 6, " or ", n.cols - 6, ").",
               call. = FALSE)
      }
      else
        return(x)
    }
    ##
    if (method == "RoB2") {
      if (n.domains != n.cols) {
        if (n.domains == n.cols - 5) {
          x <- paste0("(", LETTERS[5 + seq_len(n.cols - 5)], ") ", x)
          return(c(dm, x))
        }
        else if (n.domains < n.cols) {
          x <- paste0("(", LETTERS[5 + seq_len(n.cols - 5)], ") ",
                      "Additional item ", seq_len(n.cols - 5))
          return(c(dm, x))
        }
        else
          stop("Wrong number of domains provided for '", method,
               "' (must be ", 5, " or ", n.cols - 5, ").",
               call. = FALSE)
      }
      else
        return(x)
    }
    ##
    if (method %in% c("ROBINS-I", "ROBINS-E")) {
      if (n.domains != n.cols) {
        if (n.domains == n.cols - 7) {
          x <- paste0("(", LETTERS[7 + seq_len(n.cols - 7)], ") ", x)
          return(c(dm, x))
        }
        else if (n.domains < n.cols) {
          x <- paste0("(", LETTERS[7 + seq_len(n.cols - 7)], ") ",
                      "Additional item ", seq_len(n.cols - 7))
          return(c(dm, x))
        }
        else
          stop("Wrong number of domains provided for '", method,
               "' (must be ", 7, " or ", n.cols - 7, ").",
               call. = FALSE)
      }
      else
        return(x)
    }
    
    x
  }

  
  ##
  ##
  ## (4) Risk of bias categories
  ##
  ##
  
  if (is.null(categories)) {
    if (method == "RoB1")
      categories <-
        c("Low risk of bias", "Unclear risk of bias", "High risk of bias")
    else if (is.RoB)
      categories <-
        c("Low risk of bias", "Some concerns", "High risk of bias")
    else if (is.ROBINS)
      categories <-
        c("Low risk", "Some concerns", "High risk", "Very high risk", "NI")
    else
      stop("Argument 'categories' must be specified for unknown RoB method.",
           call. = FALSE)
  }
  else {
    if (is.RoB & length(categories) != 3)
      stop("Three categories must be provided if method = '", method, "'.",
           call. = FALSE)
    if (is.ROBINS & length(categories) != 5)
      stop("Five categories must be provided if method = '", method, "'.",
           call. = FALSE)
  }
  ##
  n.cat <- length(categories)
  
  
  ##
  ##
  ## (5) Risk of bias table
  ##
  ##
  
  rob <- data.frame(A = setcat(item1, categories))
  if (!is.null(item2))
    rob$B <- setcat(item2, categories)
  if (!is.null(item3))
    rob$C <- setcat(item3, categories)
  if (!is.null(item4))
    rob$D <- setcat(item4, categories)
  if (!is.null(item5))
    rob$E <- setcat(item5, categories)
  if (!is.null(item6))
    rob$F <- setcat(item6, categories)
  if (!is.null(item7))
    rob$G <- setcat(item7, categories)
  if (!is.null(item8))
    rob$H <- setcat(item8, categories)
  if (!is.null(item9))
    rob$I <- setcat(item9, categories)
  if (!is.null(item10))
    rob$J <- setcat(item10, categories)
  ##
  ##rob <- rob %>% mutate(across(colnames(rob), factor, categories))
  ##
  ## RoB symbols
  ##
  if (is.null(symbols)) {
    if (method == "user-defined")
      symbols <- FALSE
    else {
      if (is.RoB)
        symbols <- c("+", "?", "-")
      else if (is.ROBINS)
        symbols <- FALSE
    }
  }
  else {
    chkchar(symbols, nchar = 1)
    ##
    if (length(symbols) == 1 && is.logical(symbols)) {
      if (symbols) {
        if (method == "user-defined")
          symbols <- seq_along(categories)
        else if (is.RoB)
          symbols <- c("+", "?", "-")
        else
          symbols <- FALSE
      }
    }
    else 
      chklength(symbols, n.cat,
                text = "Wrong number of RoB symbols (argument 'symbols').")
  }
  ##
  ## RoB colours
  ##
  if (is.null(col)) {
    if (method == "user-defined")
      col <- seq_len(length(unique(unlist(rob))))
    else {
      if (is.RoB)
        col <- c("green", "yellow", "red")
      else if (is.ROBINS)
        col <-  c("green", "yellow", "red", "darkred", "darkgrey")
    }
  }
  else
    chklength(col, n.cat,
              text = "Wrong number of RoB colours (argument 'col').")
  ##
  ## RoB domain names
  ##
  if (is.null(domains)) {
    if (method == "user-defined")
      domains <- paste("Item", seq_len(ncol(rob)))
    else
      domains <- dm
  }
  ##
  domains <- setdom(domains, method, dm, ncol(rob))
  
  
  ##
  ##
  ## (6) Return risk of bias table
  ##
  ##
  
  rob <- cbind(studlab, rob)
  ##
  attr(rob, "method") <- method
  attr(rob, "domains") <- domains
  attr(rob, "symbols") <- symbols
  attr(rob, "col") <- col
  attr(rob, "legend") <- legend
  ##
  class(rob) <- c("rob", class(rob))
  ##
  if (is.null(meta.object))
    res <- rob
  else {
    res <- meta.object
    if (!is.null(res$rob) & !overwrite)
      stop("Risk of bias available ",
           "(use argument 'overwrite = TRUE' to replace RoB table",
           call. = FALSE)
    res$rob <- rob
  }
  
  
  res
}





#' @rdname rob
#' @method print rob
#' @export


print.rob <- function(x, legend = attr(x, "legend"), ...) {

  chkclass(x, "rob")
  chklogical(legend)
  
  x.prt <- x
  class(x.prt) <- "data.frame"
  print(x.prt, ...)
  
  if (legend) {
    txt.legend <-
      paste0("\nRisk of bias legend:",
             paste0("\n", attr(x, "domains"), collapse = ""),
             "\n")
    cat(txt.legend)
  }
  
  invisible(NULL)
}
