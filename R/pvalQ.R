pvalQ <- function(Q, df, lower.tail = FALSE) {
  ##
  if (length(df) == 1 & length(Q) > 1)
    df <- rep(df, length(Q))
  else if (length(df) > 1 & length(Q) == 1)
    Q <- rep(Q, length(df))
  else if (length(df) > 1 & length(Q) > 1 & length(df) != length(Q))
    stop("Length of arguments 'Q' and 'df' do not match.")
  ##
  res <- ifelse(is.na(df) | df < 1,
                NA,
                pchisq(Q, df, lower.tail = lower.tail))
  ##
  res
}
