calcH <- function(Q, df, level) {
  ##
  ## Calculate H
  ## Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58
  ##
  k <- df + 1
  ##
  if (!is.na(k)) {
    if (k > 1)
      H <- sqrt(Q / (k - 1))
    else
      H <- NA
    ##
    selogH <- ifelse(Q > k,
                     ifelse(k >= 2,
                            0.5 * (log(Q) - log(k - 1)) /
                            (sqrt(2 * Q) - sqrt(2 * k - 3)),
                            NA),
                     ifelse(k > 2,
                            sqrt(1 / (2 * (k - 2)) * (1 - 1 / (3 * (k - 2)^2))),
                            NA))
  }
  else {
    H <- NA
    selogH <- NA
  }
  ##
  tres <- ci(log(H), selogH, level)
  ##
  res <- list(TE = max(exp(tres$TE), 1),
              lower = max(exp(tres$lower), 1),
              upper = max(exp(tres$upper), 1))
  
  res
}
