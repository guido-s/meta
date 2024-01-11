smlab <- function(sm, backtransf, pscale, irscale) {
  res <- sm
  #
  if (backtransf) {
    if (sm == "ZCOR")
      res <- "COR"
    else if (is_mean(sm))
      res <- "mean"
    else if (is_prop(sm)) {
      if (pscale == 1)
        res <- "proportion"
      else
        res <- "events"
    }
    else if (is_rate(sm)) {
      if (irscale == 1)
        res <- "rate"
      else
        res <- "events"
    }
  }
  else {
    if (is_relative_effect(sm))
      res <- paste0("log", sm)
    else if (sm == "VE")
      res <- "logVR"
  }
  #
  res
}

