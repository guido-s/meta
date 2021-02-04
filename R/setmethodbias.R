setmethodbias <- function(x) {
  oldmethod <- setchar(x, .settings$meth4bias.old,
                       stop.at.error = FALSE)
  ##
  if (is.null(oldmethod))
    res <- setchar(x, .settings$meth4bias)
  else
    res <- switch(oldmethod,
                  rank = "Begg",
                  linreg = "Egger",
                  mm = "Thompson",
                  count = "Schwarzer",
                  score = "Harbord")
  ##
  res
}
