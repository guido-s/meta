blup_bounds_fixture <- function(sm = "VE", subgroup = FALSE, backtransf = TRUE) {
  metagen(c(-0.4, -0.1, 0.2, 0.7, -0.2, 0.3, 0.6, 0.9), rep(0.2, 8),
          sm = sm, subgroup = if (subgroup) rep(c("A", "B"), each = 4) else NULL,
          method.tau = "DL", method.random.ci = "classic", prediction = TRUE,
          backtransf = backtransf)
}
