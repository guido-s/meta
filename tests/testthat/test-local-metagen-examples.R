test_that("metagen rejects n.c for single-arm summary measures", {
  TE <- c(0.1, 0.2, 0.3)
  seTE <- c(0.01, 0.02, 0.03)
  n <- c(10, 20, 30)

  for (sm in c("PRAW", "PLOGIT", "PLN", "PAS", "PFT",
               "IR", "IRLN", "IRS", "IRFT",
               "MRAW", "MLN",
               "COR", "ZCOR")) {
    expect_error(
      metagen(TE, seTE, sm = sm, n.c = n),
      "Argument 'n.c' must not be provided"
    )
  }
})
