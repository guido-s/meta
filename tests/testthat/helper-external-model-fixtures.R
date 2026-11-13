external_multilevel_fixture <- function() {
  data.frame(
    district = rep(seq_len(6), each = 3),
    idx = seq_len(18),
    group = factor(rep(c("A", "B"), each = 9)),
    yi = rep(c(-0.6, 0.1, 0.7, -0.2, 0.5, 1.1), each = 3) +
      c(-0.35, 0.12, 0.4, 0.3, -0.4, 0.1, -0.2, 0.4, -0.1,
        0.2, -0.3, 0.35, -0.4, 0.15, 0.3, 0.35, -0.1, -0.3),
    vi = rep(c(0.025, 0.04, 0.06), 6))
}
