percent_label <- function(x, type = "CI", digits = 1,
                          func = round, lab.NA = NULL, ...) {
  value <- func(100 * x, digits = digits, ...)
  #
  if (length(type) > 0 && type != "")
    suffix <- paste0(" ", type)
  else
    suffix <- ""
  
  res <- paste0(value, "%", suffix)
  #
  if (!is.null(lab.NA))
    res[is.na(x)] <- lab.NA
  #
  res
}
