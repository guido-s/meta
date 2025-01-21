#' Produce traffic light plot of risk of bias assessment
#' 
#' @description
#' Produce traffic light plot of risk of bias assessment
#' 
#' @param object An object of class \code{rob}.
#' @param colour Specify colour scheme for the traffic light plot; see
#'   \code{\link[robvis]{rob_summary}}.
#' @param psize Size of the traffic lights.
#' @param quiet A logical to suppress the display of the traffic light
#'   plot.
#' 
#' @details
#' This is a wrapper function for
#' \code{\link[robvis]{rob_traffic_light}} of R package \bold{robvis}
#' to produce a traffic light plot of risk of bias assessment.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{rob}}, \code{\link{barplot.rob}},
#'   \code{\link[robvis]{rob_traffic_light}}
#' 
#' @keywords hplot
#' 
#' @examples
#' # Use RevMan 5 settings
#' oldset <- settings.meta("RevMan5")
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
#' \dontrun{
#' # Traffic light plot (R package 'robvis' must be available)
#' if (requireNamespace("robvis", quietly = TRUE))
#'  traffic_light(rob(m2))
#' }
#' 
#' # Use previous settings
#' settings.meta(oldset)
#' 
#' @export traffic_light


traffic_light <- function(object,
                          colour = "cochrane",
                          psize = 15,
                          quiet = FALSE) {
  
  chkclass(object, "rob")
  rob <- object
  ##
  tool <- attr(object, "tool")
  ##
  chklogical(quiet)
  
  
  if (!is_installed_package("robvis", stop = FALSE))
    stop(paste0("Package 'robvis' missing.",
                "\n  ",
                "Please use the following R command for installation:",
                "\n  install.packages(\"robvis\")"),
         call. = FALSE)
  
  
  if (!(tool %in% c("RoB1", "RoB2", "ROBINS-I", "user-defined"))) {
    warning("R function 'rob_traffic_light() in R package robvis not usable ",
            "with risk of bias tool '", tool, "'.",
            call. = FALSE)
    return(invisible(NULL))
  }
  else {
    if (tool %in% c("RoB2", "ROBINS-I") & is.null(rob$Overall)) {
      warning("No traffic light plot created as ",
              "overall risk of bias assessment is missing.",
              call. = FALSE)
      return(invisible(NULL))
    }
    ##
    if (is.null(rob$Weight))
      rob$Weight <- 1
    ##
    if (tool %in% c("RoB1", "user-defined")) {
      domains <- attr(rob, "domains")
      ##
      nam <- names(rob)
      ##
      if (tool == "RoB1") {
        nam[nam == "A"] <- "Random.sequence.generation."
        nam[nam == "B"] <- "Allocation.concealment."
        nam[nam == "C"] <- "Blinding.of.participants.and.personnel."
        nam[nam == "D"] <- "Blinding.of.outcome.assessment."
        nam[nam == "E"] <- "Incomplete.outcome.data."
        nam[nam == "F"] <- "Selective.reporting."
        nam[nam == "G"] <- "Other.bias."
        ##
        if (length(domains) > 8)
          nam[nam == "H"] <- paste0(domains[8], ".")
        if (length(domains) > 9)
          nam[nam == "I"] <- paste0(domains[9], ".")
        if (length(domains) > 10)
          nam[nam == "J"] <- paste0(domains[10], ".")
      }
      else if (tool == "user-defined") {
        domains[domains == "Overall risk of bias"] <- "Overall"
        ##
        nam[nam == "Overall"] <- "O"
        nam[nam %in% LETTERS[1:15]] <- paste0(gsub(" ", ".", domains), ".")
        ##
        tool <- "RoB1"
      }
      ##
      names(rob) <- nam
    }
    ##
    return(robvis::rob_traffic_light(rob, tool = toupper(tool),
                                     colour = colour, psize = psize,
                                     quiet = quiet))
  }
}
