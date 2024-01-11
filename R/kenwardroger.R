kenwardroger <- function(w.random) {
  #
  # Kenward-Roger method for meta-analysis
  # (Partlett & Riley, 2017, Stat Med)
  #
  w1p <- sum(w.random)
  w2p <- sum(w.random^2)
  w3p <- sum(w.random^3)
  #
  IE <- w2p * 0.5 - w3p / w1p + 0.5 * (w2p / w1p)^2
  var <- (1 + 2 * (w3p / w1p - (w2p / w1p)^2) / IE) / w1p
  df <- 2 * IE / (var * w2p)^2
  #
  res <- list(se = if (!is.na(var) && var >= 0) sqrt(var) else NA, df = df,
              IE = IE, w1p = w1p, w2p = w2p, w3p = w3p)
  res
}
