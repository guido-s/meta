test_that("sequential prediction defaults use DL random-effects results", {
  TE <- c(-0.8, 0.9, -0.4, 1.2, 0.3, 1.6)
  seTE <- c(0.12, 0.3, 0.18, 0.35, 0.22, 0.4)
  fields <- c("TE", "seTE", "lower", "upper", "tau2",
              "lower.predict", "upper.predict", "df.predict")

  for (algorithm in c("metacum", "metainf")) {
    sequential <- get(algorithm)
    for (inherited in c(FALSE, TRUE)) {
      ma <- metagen(TE, seTE, common = TRUE, random = FALSE,
                    prediction = inherited, method.tau = "DL",
                    method.random.ci = "classic", method.predict = "V",
                    level.ma = 0.95, level.predict = 0.95)
      original <- serialize(ma, NULL)
      selected <- if (inherited) sequential(ma) else
        sequential(ma, prediction = TRUE)
      explicit <- sequential(ma, pooled = "random", prediction = TRUE)
      common <- sequential(ma, prediction = FALSE)
      for (field in c(fields, paste0(fields[1:4], ".pooled"),
                      "lower.predict.pooled", "upper.predict.pooled"))
        expect_equal(selected[[field]], explicit[[field]])

      # Check every prefix / leave-one-out subset, and the original pooled fit.
      for (i in seq_len(length(TE) + 1L)) {
        pooled <- i > length(TE)
        sel <- if (pooled) seq_along(TE) else if (algorithm == "metacum")
          seq_len(i) else seq_along(TE)[-i]
        n <- length(sel)
        y <- TE[sel]
        v <- seTE[sel]^2
        w <- 1 / v
        mu.common <- sum(w * y) / sum(w)
        tau2 <- if (n == 1L) 0 else
          max(0, (sum(w * (y - mu.common)^2) - (n - 1)) /
                (sum(w) - sum(w^2) / sum(w)))
        wr <- 1 / (v + tau2)
        mu <- sum(wr * y) / sum(wr)
        se <- sqrt(1 / sum(wr))
        expected <- c(TE = mu, seTE = se,
                      lower = mu - qnorm(0.975) * se,
                      upper = mu + qnorm(0.975) * se)
        for (field in names(expected)) {
          actual <- if (pooled) selected[[paste0(field, ".pooled")]] else
            selected[[field]][i]
          expect_equal(actual, unname(expected[field]), tolerance = 1e-10)
        }
        common.TE <- if (pooled) common$TE.pooled else common$TE[i]
        expect_equal(common.TE, mu.common, tolerance = 1e-10)

        # PIs are unavailable for the smallest subsets; check them from n = 3.
        if (n >= 3L) {
          expect_gt(tau2, 0)
          expect_gt(abs(mu - mu.common), 0.01)
          halfwidth <- qt(0.975, df = n - 1) * sqrt(tau2 + se^2)
          for (field in c("lower.predict", "upper.predict")) {
            actual <- if (pooled) selected[[paste0(field, ".pooled")]] else
              selected[[field]][i]
            expected.limit <- mu + if (field == "lower.predict")
              -halfwidth else halfwidth
            expect_equal(actual, expected.limit, tolerance = 1e-10)
          }
        }
      }
      expect_identical(serialize(ma, NULL), original)
    }
  }
})
