#' Risk of bias assessment
#' 
#' @description
#' Create table with risk of bias assessment or add table to existing
#' meta-analysis
#' 
#' @param item1 Risk of bias item 1 or a meta-analysis object of class
#'   \code{meta} with information on risk of bias assessment.
#' @param item2 Risk of bias item 2.
#' @param item3 Risk of bias item 3.
#' @param item4 Risk of bias item 4.
#' @param item5 Risk of bias item 5.
#' @param item6 Risk of bias item 6.
#' @param item7 Risk of bias item 7.
#' @param item8 Risk of bias item 8.
#' @param item9 Risk of bias item 9.
#' @param item10 Risk of bias item 10.
#' @param overall Overall risk of bias assess.
#' @param weight Weight for each study.
#' @param studlab Study labels.
#' @param data A data frame or a meta-analysis object of class
#'   \code{meta}.
#' @param tool Risk of bias (RoB) tool.
#' @param domains A character vector with names of RoB domains.
#' @param categories Possible RoB categories.
#' @param cat1 Possible categories for RoB item 1.
#' @param cat2 Possible categories for RoB item 2.
#' @param cat3 Possible categories for RoB item 3.
#' @param cat4 Possible categories for RoB item 4.
#' @param cat5 Possible categories for RoB item 5.
#' @param cat6 Possible categories for RoB item 6.
#' @param cat7 Possible categories for RoB item 7.
#' @param cat8 Possible categories for RoB item 8.
#' @param cat9 Possible categories for RoB item 9.
#' @param cat10 Possible categories for RoB item 10.
#' @param cat.overall Possible categories for overall RoB.
#' @param col Colours for RoB categories.
#' @param col1 Colours for categories for RoB item 1.
#' @param col2 Colours for categories for RoB item 2.
#' @param col3 Colours for categories for RoB item 3.
#' @param col4 Colours for categories for RoB item 4.
#' @param col5 Colours for categories for RoB item 5.
#' @param col6 Colours for categories for RoB item 6.
#' @param col7 Colours for categories for RoB item 7.
#' @param col8 Colours for categories for RoB item 8.
#' @param col9 Colours for categories for RoB item 9.
#' @param col10 Colours for categories for RoB item 10.
#' @param col.overall Colours for categories for overall RoB.
#' @param symbols Corresponding symbols for RoB categories.
#' @param symb1 Corresponding symbols for RoB item 1.
#' @param symb2 Corresponding symbols for RoB item 2.
#' @param symb3 Corresponding symbols for RoB item 3.
#' @param symb4 Corresponding symbols for RoB item 4.
#' @param symb5 Corresponding symbols for RoB item 5.
#' @param symb6 Corresponding symbols for RoB item 6.
#' @param symb7 Corresponding symbols for RoB item 7.
#' @param symb8 Corresponding symbols for RoB item 8.
#' @param symb9 Corresponding symbols for RoB item 9.
#' @param symb10 Corresponding symbols for RoB item 10.
#' @param symb.overall Corresponding symbols for overall RoB.
#' @param legend A logical specifying whether legend with RoB domains
#'   should be printed.
#' @param overwrite A logical indicating whether an existing risk of
#'   bias table in a meta-analysis object should be overwritten.
#' @param warn A logical indicating whether warnings should be
#'   printed.
#' @param x An object of class \code{rob}.
#' @param details A logical indicating whether to print details on
#'   categories and colours.
#' @param \dots Additional printing arguments.
#' 
#' @details
#' This function can be used to define a risk of bias (RoB) assessment
#' for a meta-analysis which can be shown in a forest plot
#' (\code{\link{forest.meta}}), summary weighted barplot
#' (\code{\link{barplot.rob}}) or traffic light plot
#' (\code{\link{traffic_light}}). It is also possible to extract the
#' risk of bias assessment from a meta-analysis with RoB information.
#'
#' The risk of bias table contains
#'
#' \itemize{
#' \item study labels;
#' \item variables for individual RoB domains (with variable names A,
#'   B, \dots);
#' \item an overall RoB assessment if argument \code{overall} is
#'   provided;
#' \item weights for individual studies used in summary weighted
#'   barplots.
#' }
#' Note, an overall RoB assessment is mandatory to create a summary
#' weighted barplot or a traffic light plot.
#'
#' The RoB table is directly returned if argument \code{data} is a
#' data frame or argument \code{item1} is a meta-analysis with risk of
#' bias assessment. The RoB table is added as a new list element 'rob'
#' to a meta-analysis object if argument \code{data} is a
#' meta-analysis.
#'
#' The user must either specify the categories and (optionally)
#' domains of the RoB tool (using the eponymous arguments) or one of
#' the following RoB tools.
#' \tabular{ll}{
#' \bold{Argument} \tab \bold{Risk of bias tool} \cr
#' \code{tool = "RoB1"} \tab RoB 1 tool for randomized studies
#'   (Higgins et al., 2011) \cr
#' \code{tool = "RoB2"} \tab RoB 2 tool for randomized studies
#'   (Sterne et al., 2019) \cr
#' \code{tool = "RoB2-cluster"} \tab RoB 2 tool for
#'   cluster-randomized trials \cr
#' \code{tool = "RoB2-crossover"} \tab RoB 2 tool for crossover
#'   trials \cr
#' \code{tool = "ROBINS-I"} \tab Risk Of Bias In Non-randomized
#'   Studies - of Interventions \cr
#' \tab (Sterne et al., 2016) \cr
#' \code{tool = "ROBINS-E"} \tab Risk Of Bias In Non-randomized
#'   Studies - of Exposures \cr
#' \tab (ROBINS-E Development Group, 2023)
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
#'  \item Blinding of outcome assessment (detection bias)
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
#'  \item First item
#'  \item Second item
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
#' they must be of the same length as the number of categories. 
#' }
#' 
#' @return
#' A data frame with study labels and risk of bias items and
#' additional class "rob".
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.meta}}, \code{\link{barplot.rob}},
#'   \code{\link{traffic_light}}
#'
#' @references
#' Higgins JPT, Altman DG, Gøtzsche PC, Jüni P, Moher D, Oxman AD et
#' al. (2011):
#' The Cochrane Collaboration's tool for assessing risk of bias in
#' randomised trials.
#' \emph{British Medical Journal}, \bold{343}: d5928
#'
#' ROBINS-E Development Group (Higgins J, Morgan R, Rooney A et al.) (2023):
#' Risk Of Bias In Non-randomized Studies - of Exposure (ROBINS-E)
#' Available from: \url{https://www.riskofbias.info/welcome/robins-e-tool}.
#' 
#' Sterne JA, Hernán MA, Reeves BC, Savović J, Berkman ND, Viswanathan
#' M, et al. (2016):
#' ROBINS-I: a tool for assessing risk of bias in non-randomised
#' studies of interventions.
#' \emph{British Medical Journal}, \bold{355}: i4919
#' 
#' Sterne JAC, Savović J, Page MJ, Elbers RG, Blencowe NS, Boutron I,
#' et al. (2019):
#' RoB 2: a revised tool for assessing risk of bias in randomised trials.
#' \emph{British Medical Journal}, \bold{366}: l4898. 
#' 
#' @examples
#' # Use RevMan 5 settings
#' oldset <- settings.meta("RevMan5", quietly = FALSE)
#' 
#' data(caffeine)
#' 
#' m1 <- metabin(h.caf, n.caf, h.decaf, n.decaf, sm = "OR",
#'   data = caffeine, studlab = paste(study, year))
#'
#' # Add risk of bias assessment to meta-analysis
#' m2 <- rob(D1, D2, D3, D4, D5, overall = rob, data = m1, tool = "rob2")
#' 
#' # Print risk of bias assessment
#' rob(m2)
#' 
#' # Forest plot with risk of bias assessment
#' forest(m2)
#'
#' # Use previous settings
#' settings.meta(oldset)
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
                overall = NULL,
                weight = NULL,
                data = NULL,
                ##
                tool = gs("tool.rob"),
                domains = NULL,
                ##
                categories = NULL,
                cat1 = categories,
                cat2 = categories,
                cat3 = categories,
                cat4 = categories,
                cat5 = categories,
                cat6 = categories,
                cat7 = categories,
                cat8 = categories,
                cat9 = categories,
                cat10 = categories,
                cat.overall = categories,
                ##
                col = NULL,
                col1 = col,
                col2 = col,
                col3 = col,
                col4 = col,
                col5 = col,
                col6 = col,
                col7 = col,
                col8 = col,
                col9 = col,
                col10 = col,
                col.overall = col,
                ##
                symbols = NULL,
                symb1 = symbols,
                symb2 = symbols,
                symb3 = symbols,
                symb4 = symbols,
                symb5 = symbols,
                symb6 = symbols,
                symb7 = symbols,
                symb8 = symbols,
                symb9 = symbols,
                symb10 = symbols,
                symb.overall = symbols,
                ##
                legend = TRUE,
                ##
                overwrite = FALSE,
                warn = TRUE) {
  
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
    ##
    data <- data$data
    if (!is.null(meta.object$subset))
      data <- data[meta.object$subset, ]
    ##
    if (is.null(studlab))
      studlab <- data$.studlab    
  }
  ##
  ## Catch 'item1', etc.
  ##
  item1 <- catch("item1", mc, data, sfsp)
  if (!is.null(item1) && inherits(item1, "meta")) {
    if (is.null(item1$rob))
      return(NULL)
    else
      return(item1$rob)
  }
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
  missing.overall <- missing(overall)
  overall <- catch("overall", mc, data, sfsp)
  avail.overall <- !missing.overall & !is.null(overall)
  ##
  missing.weight <- missing(weight)
  weight <- catch("weight", mc, data, sfsp)
  avail.weight <- !missing.weight & !is.null(weight)
  ##
  studlab <- catch("studlab", mc, data, sfsp)
  
  
  ##
  ##
  ## (2) Check data
  ##
  ##
  
  fun <- "rob"
  ##
  avail1 <- !is.null(item1)
  avail2 <- !is.null(item2)
  avail3 <- !is.null(item3)
  avail4 <- !is.null(item4)
  avail5 <- !is.null(item5)
  avail6 <- !is.null(item6)
  avail7 <- !is.null(item7)
  avail8 <- !is.null(item8)
  avail9 <- !is.null(item9)
  avail10 <- !is.null(item10)
  ##
  if (avail1)
    k.All <- length(item1)
  else if (avail2)
    k.All <- length(item2)
  else if (avail3)
    k.All <- length(item3)
  else if (avail4)
    k.All <- length(item4)
  else if (avail5)
    k.All <- length(item5)
  else if (avail6)
    k.All <- length(item6)
  else if (avail7)
    k.All <- length(item7)
  else if (avail8)
    k.All <- length(item8)
  else if (avail9)
    k.All <- length(item9)
  else if (avail10)
    k.All <- length(item10)
  else
    stop("No information on risk of bias domains provided.", call. = FALSE)
  ##
  if (avail2)
    chklength(item2, k.All, fun)
  ##
  if (avail3)
    chklength(item3, k.All, fun)
  ##
  if (avail4)
    chklength(item4, k.All, fun)
  ##
  if (avail5)
    chklength(item5, k.All, fun)
  ##
  if (avail6)
    chklength(item6, k.All, fun)
  ##
  if (avail7)
    chklength(item7, k.All, fun)
  ##
  if (avail8)
    chklength(item8, k.All, fun)
  ##
  if (avail9)
    chklength(item9, k.All, fun)
  ##
  if (avail10)
    chklength(item10, k.All, fun)
  ##
  if (avail.overall)
    chklength(overall, k.All, fun)
  ##
  if (avail.weight & length(weight) > 1)
    chklength(weight, k.All, fun)
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
  chklogical(warn)
  ##
  if (!is.null(tool))
    tool <- setchar(tool, gs("tool4rob"))
  else
    tool <- "user-defined"
  ##
  
  
  ##
  ##
  ## (3) Risk of bias categories
  ##
  ##
  cat1 <- definecat(avail1, cat1, item1, tool, warn)
  cat2 <- definecat(avail2, cat2, item2, tool, warn)
  cat3 <- definecat(avail3, cat3, item3, tool, warn)
  cat4 <- definecat(avail4, cat4, item4, tool, warn)
  cat5 <- definecat(avail5, cat5, item5, tool, warn)
  cat6 <- definecat(avail6, cat6, item6, tool, warn)
  cat7 <- definecat(avail7, cat7, item7, tool, warn)
  cat8 <- definecat(avail8, cat8, item8, tool, warn)
  cat9 <- definecat(avail9, cat9, item9, tool, warn)
  cat10 <- definecat(avail10, cat10, item10, tool, warn)
  ##
  cat.overall <- definecat(avail.overall, cat.overall, overall, tool, warn)
  
  
  ##
  ##
  ## (4) Create risk of bias table
  ##
  ##
  
  rob <- data.frame(Study = studlab)
  ##
  if (avail1)
    rob$A <- setcat(item1, cat1)
  ##
  if (avail2)
    rob$B <- setcat(item2, cat2)
  ##
  if (avail3)
    rob$C <- setcat(item3, cat3)
  ##
  if (avail4)
    rob$D <- setcat(item4, cat4)
  ##
  if (avail5)
    rob$E <- setcat(item5, cat5)
  ##
  if (avail6)
    rob$F <- setcat(item6, cat6)
  ##
  if (avail7)
    rob$G <- setcat(item7, cat7)
  ##
  if (avail8)
    rob$H <- setcat(item8, cat8)
  ##
  if (avail9)
    rob$I <- setcat(item9, cat9)
  ##
  if (avail10)
    rob$J <- setcat(item10, cat10)
  ##
  domain.available <- c(avail1, avail2, avail3, avail4, avail5, avail6, avail7,
                        avail8, avail9, avail10)
  names(domain.available) <- paste0("D", 1:10)
  ##
  if (avail.overall)
    rob$Overall <- setcat(overall, cat.overall)
  ##
  if (avail.weight && length(weight) == 1 && is.character(weight)) {
    if (is.meta)
      weight <- setchar(weight, c("common", "random"))
    else
      weight <- NULL
  }
  ##
  if (!avail.weight && is.meta && meta.object$overall)
    weight <-
      if (meta.object$random)
        "random"
      else if (meta.object$common)
        "common"
      else
        NULL
  ##
  if (!is.null(weight)) {
    if (length(weight) == 1 && is.character(weight)) {
      sel.w <-
        if (is.null(meta.object$subset))
          seq_along(meta.object$TE)
        else
          meta.object$subset
      ##
      if (weight == "random")
        rob$Weight <- meta.object$w.random[sel.w]
      else if (weight == "common")
        rob$Weight <- meta.object$w.common[sel.w]
    }
    else
      rob$Weight <- weight
    ##
    avail.weight <- TRUE
  }
  else
    avail.weight <- FALSE
  
  
  ##
  ##
  ## (4) Names for risk of bias domains
  ##
  ##
  
  if (tool == "user-defined")
    dm <- paste("Item", 1:10)
  else if (tool == "RoB1")
    dm <-
      c("Random sequence generation (selection bias)",
        "Allocation concealment (selection bias)",
        "Blinding of participants and personnel (performance bias)",
        "Blinding of outcome assessment (detection bias)",
        "Incomplete outcome data (attrition bias)",
        "Selective reporting (reporting bias)",
        "Other bias",
        rep("", 10 - 7))
  else if (tool == "RoB2")
    dm <-
      c("Bias arising from the randomization process",
        "Bias due to deviations from intended intervention",
        "Bias due to missing outcome data",
        "Bias in measurement of the outcome",
        "Bias in selection of the reported result",
        rep("", 10 - 5))
  else if (tool == "RoB2-cluster")
    dm <-
      c("Bias arising from the randomization process",
        paste("Bias arising from the identification or",
              "recruitment of participants into clusters"),
        "Bias due to deviations from intended intervention",
        "Bias due to missing outcome data",
        "Bias in measurement of the outcome",
        "Bias in selection of the reported result",
        rep("", 10 - 6))
  else if (tool == "RoB2-crossover")
    dm <-
      c("Bias arising from the randomization process",
        "Bias arising from period and carryover effects",
        "Bias due to deviations from intended intervention",
        "Bias due to missing outcome data",
        "Bias in measurement of the outcome",
        "Bias in selection of the reported result",
        rep("", 10 - 6))
  else if (tool == "ROBINS-I")
    dm <-
      c("Risk of bias due to confounding",
        "Risk of bias in selection of participants into the study",
        "Risk of bias in classification of interventions",
        "Risk of bias due to deviations from intented interventions",
        "Risk of bias due to missing outcome data",
        "Risk of bias in measurement of the outcome",
        "Risk of bias in the selection of the reported results",
        rep("", 10 - 7))
  else if (tool == "ROBINS-E")
    dm <-
      c("Risk of bias due to confounding",
        "Risk of bias arising from measurement of the exposure",
        paste("Risk of bias in selection of participants",
              "into the study (or into the analysis)"),
        "Risk of bias due to post-exposure interventions",
        "Risk of bias due to missing data",
        "Risk of bias in measurement of the outcome",
        "Risk of bias in selection of the reported results",
        rep("", 10 - 7))
  else
    dm <- rep("", 10)
  ##
  if (tool %in% c("RoB1", "ROBINS-I", "ROBINS-E")) {
    if (avail8)
      dm[8] <- "Additional item 1"
    ##
    if (avail9)
      dm[9] <- "Additional item 2"
    ##
    if (avail10)
      dm[10] <- "Additional item 3"
  }
  else if (tool %in% c("RoB2-cluster", "RoB2-crossover")) {
    if (avail7)
      dm[7] <- "Additional item 1"
    ##
    if (avail8)
      dm[8] <- "Additional item 2"
    ##
    if (avail9)
      dm[9] <- "Additional item 3"
    ##
    if (avail10)
      dm[10] <- "Additional item 4"
  }
  else if (tool == "RoB2") {
    if (avail6)
      dm[6] <- "Additional item 1"
    ##
    if (avail7)
      dm[7] <- "Additional item 2"
    ##
    if (avail8)
      dm[8] <- "Additional item 3"
    ##
    if (avail9)
      dm[9] <- "Additional item 4"
    ##
    if (avail10)
      dm[10] <- "Additional item 5"
  }
  ##
  if (is.null(domains))
    domains <- dm[domain.available]
  ##
  else if (tool == "user-defined") {
    if (length(domains) != sum(domain.available))
      stop("Number of domain names does not match number of domains.",
           call. = FALSE)
  }
  else
    domains <- setdom(dm, tool, domains, domain.available)
  ##
  if (avail.overall)
    domains <- c(domains, "Overall risk of bias")
  
  
  ##
  ##
  ## (6) Risk of bias symbols
  ##
  ##
  
  symb1 <- definesymb(avail1, symb1, cat1, tool, warn)
  symb2 <- definesymb(avail2, symb2, cat2, tool, warn)
  symb3 <- definesymb(avail3, symb3, cat3, tool, warn)
  symb4 <- definesymb(avail4, symb4, cat4, tool, warn)
  symb5 <- definesymb(avail5, symb5, cat5, tool, warn)
  symb6 <- definesymb(avail6, symb6, cat6, tool, warn)
  symb7 <- definesymb(avail7, symb7, cat7, tool, warn)
  symb8 <- definesymb(avail8, symb8, cat8, tool, warn)
  symb9 <- definesymb(avail9, symb9, cat9, tool, warn)
  symb10 <- definesymb(avail10, symb10, cat10, tool, warn)
  ##
  symb.overall <-
    definesymb(avail.overall, symb.overall, cat.overall, tool, warn)
  
  
  ##
  ##
  ## (7) Risk of bias colours
  ##
  ##
  
  col1 <- definecol(avail1, col1, cat1, tool, warn)
  col2 <- definecol(avail2, col2, cat2, tool, warn)
  col3 <- definecol(avail3, col3, cat3, tool, warn)
  col4 <- definecol(avail4, col4, cat4, tool, warn)
  col5 <- definecol(avail5, col5, cat5, tool, warn)
  col6 <- definecol(avail6, col6, cat6, tool, warn)
  col7 <- definecol(avail7, col7, cat7, tool, warn)
  col8 <- definecol(avail8, col8, cat8, tool, warn)
  col9 <- definecol(avail9, col9, cat9, tool, warn)
  col10 <- definecol(avail10, col10, cat10, tool, warn)
  ##
  col.overall <- definecol(avail.overall, col.overall, cat.overall, tool, warn)
  
  
  ##
  ##
  ## (8) Return risk of bias table
  ##
  ##
  
  attr(rob, "tool") <- tool
  attr(rob, "domains") <- domains
  ##
  attr(rob, "overall") <- avail.overall
  attr(rob, "weight") <- avail.weight
  ##
  attr(rob, "categories") <-
    compact(list(cat1 = cat1, cat2 = cat2, cat3 = cat3, cat4 = cat4,
                 cat5 = cat5, cat6 = cat6, cat7 = cat7, cat8 = cat8,
                 cat9 = cat9, cat10 = cat10, cat.overall = cat.overall))
  attr(rob, "symbols") <-
    compact(list(symb1 = symb1, symb2 = symb2, symb3 = symb3, symb4 = symb4,
                 symb5 = symb5, symb6 = symb6, symb7 = symb7, symb8 = symb8,
                 symb9 = symb9, symb10 = symb10, symb.overall = symb.overall))
  attr(rob, "col") <-
    compact(list(col1 = col1, col2 = col2, col3 = col3, col4 = col4,
                 col5 = col5, col6 = col6, col7 = col7, col8 = col8,
                 col9 = col9, col10 = col10, col.overall = col.overall))
  ##
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


print.rob <- function(x, legend = attr(x, "legend"), details = TRUE, ...) {

  chkclass(x, "rob")
  ##
  chklogical(legend)
  chklogical(details)
  
  x.prt <- x
  class(x.prt) <- "data.frame"
  
  print(x.prt, ...)
  
  if (legend) {
    cat(
      paste0("\nRisk of bias domains:",
             paste0("\n", catleg(x), collapse = ""),
             "\n")
    )
    ##
    if (details) {
      cat(
        paste0("\nCategories:",
               paste0("\n", catcat(x), collapse = ""),
               "\n")
      )
      ##
      cat(
        paste0("\nColours:",
               paste0("\n", catcol(x), collapse = ""),
               "\n")
      )
      ##
      symbols <- unique(unlist(attr(x, "symbols")))
      if (!(length(symbols) == 1 && is.logical(symbols) && !symbols))
      cat(
        paste0("\nSymbols:",
               paste0("\n", catsymb(x), collapse = ""),
               "\n")
      )        
    }
  }
  
  
  invisible(NULL)
}
