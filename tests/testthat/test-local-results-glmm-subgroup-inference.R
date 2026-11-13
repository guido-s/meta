test_that("GLMM subgroup coefficients use joint-fit DF in either method order", {
  singles <- lapply(c("classic", "HK"), subgroup_inference_fixture)
  names(singles) <- c("classic", "HK")
  fits <- c(singles, list(subgroup_inference_fixture(c("classic", "HK")),
                         subgroup_inference_fixture(c("HK", "classic"))))
  fields <- c("seTE.random.w", "df.random.w", "statistic.random.w",
              "pval.random.w", "lower.random.w", "upper.random.w")
  for (x in fits) {
    for (method in x$method.random.ci) {
      column <- match(method, x$method.random.ci)
      values <- lapply(x[fields], function(value) {
        if (is.matrix(value)) unname(value[, column]) else unname(value)
      })
      reference <- lapply(singles[[method]][fields], unname)
      expect_equal(values, reference)
      expect_equal(x$TE.random.w, singles[[method]]$TE.random.w)

      # Both coefficients come from one ten-study, two-coefficient fit.
      df <- if (method == "HK") length(x$TE) - 2 else Inf
      expect_equal(values$df.random.w, rep(df, 2))
      statistic <- unname(x$TE.random.w) / values$seTE.random.w
      expect_equal(values$statistic.random.w, statistic)
      expect_equal(values$pval.random.w,
                   2 * pt(abs(statistic), df = values$df.random.w,
                          lower.tail = FALSE))
      margin <- qt((1 + x$level.ma) / 2, df = values$df.random.w) *
        values$seTE.random.w
      expect_equal(values$lower.random.w, unname(x$TE.random.w) - margin)
      expect_equal(values$upper.random.w, unname(x$TE.random.w) + margin)
    }
  }
})
