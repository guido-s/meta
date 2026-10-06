smlab <- function(sm, backtransf, pscale, irscale,
                  func.backtransf = NULL, forest = FALSE) {
  
  res <- sm
  #
  if (backtransf) {
    if (is_cor(sm))
      res <- "correlation"
    #
    else if (is_mean(sm))
      res <- "mean"
    #
    else if (is_prop(sm) || sm == "proportion") {
      if (pscale == 1)
        res <- "proportion"
      else
        res <- "events"
    }
    #
    else if (is_rate(sm)) {
      if (irscale == 1)
        res <- "rate"
      else
        res <- "events"
    }
    #
    # Print first letter in upper case
    #
    res <- str_replace(res, "^.", str_to_upper)
  }
  else {
    if (sm == "VE")
      res <- paste0(gs("log.prefix"), "VR")
    else if (is_relative_effect(sm) ||
             (!is.null(func.backtransf) && func.backtransf == "exp"))
      res <- paste0(gs("log.prefix"), sm)
  }
  #
  res
}
