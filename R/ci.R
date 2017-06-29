ci <- function(TE, seTE, level = 0.95, df = NULL, null.effect = 0) {

  if (level <= 0 | level >= 1)
    stop("no valid level for confidence interval")

  alpha <- 1 - level
  
  if (is.null(df)) {
    lower  <- TE - qnorm(1 - alpha / 2)*seTE
    upper  <- TE + qnorm(1 - alpha / 2)*seTE
    zscore <- (TE - null.effect) / seTE
    pval   <- 2 * pnorm(abs(zscore), lower.tail = FALSE)
    df <- NA
  }
  else {
    df[df == 0] <- NA
    ##
    zscore <- (TE - null.effect) / seTE
    ##
    if (length(df) == 1) {
      if (is.na(df)) {
        lower <- TE - qnorm(1 - alpha / 2) * seTE
        upper <- TE + qnorm(1 - alpha / 2) * seTE
        ##
        pval <- 2 * pnorm(abs(zscore), lower.tail = FALSE)
      }
      else {
        lower <- TE - qt(1 - alpha / 2, df = df) * seTE
        upper <- TE + qt(1 - alpha / 2, df = df) * seTE
        ##
        pval <- 2 * pt(abs(zscore), df = df, lower.tail = FALSE)
      }
    }
    else {
      lower <- ifelse(!is.na(df),
                      TE - qt(1 - alpha / 2, df = df) * seTE,
                      TE - qnorm(1 - alpha / 2) * seTE)
      upper <- ifelse(!is.na(df),
                      TE + qt(1 - alpha / 2, df = df) * seTE,
                      TE + qnorm(1 - alpha / 2) * seTE)
      ##        
      pval <- ifelse(!is.na(df),
                     2 * pt(abs(zscore), df = df, lower.tail = FALSE),
                     2 * pnorm(abs(zscore), lower.tail = FALSE))
    }
    ##
    df[is.na(df)] <- 0
  }

  list(TE = TE, seTE = seTE,
       lower = lower, upper = upper,
       z = zscore, p = pval, level = level,
       df = df, null.effect = null.effect)
}
