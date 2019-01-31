seTE.pval <- function(TE, pval, df = rep_len(NA, length(TE))) {
  
  
  ##
  ## Check arguments
  ##
  if (missing(TE))
    stop("Mandatory argument 'TE' missing.", call. = FALSE)
  if (missing(pval))
    stop("Mandatory argument 'pval' missing", call. = FALSE)
  ##
  k <- length(TE)
  arg <- "TE"
  chklength(pval, k, arg)
  chklength(df, k, arg)
  ##
  if (any(pval <= 0, na.rm = TRUE) | any(pval >= 1, na.rm = TRUE))
    stop("No valid value for p-value", call. = FALSE)
  
  
  ##
  ## Parmar et al. (1998), Stat Med
  ##
  ## Equation (7)
  ##
  varTE <- ifelse(is.na(df),
                  (TE / qnorm(pval / 2, lower.tail = FALSE))^2,
                  (TE / qt(pval / 2, df = df, lower.tail = FALSE))^2)
  seTE <- sqrt(varTE)
  
  
  res <- list(TE = TE, seTE = seTE, pval = pval)
  ##
  res
}
