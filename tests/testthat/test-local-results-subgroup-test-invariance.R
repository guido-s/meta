subgroup_test_invariance_fixture <- function(methods, level = 0.95,
                                               tau.common = FALSE, adhoc = "",
                                               low.heterogeneity = FALSE) {
  TE <- if (low.heterogeneity)
    c(0.1, 0.11, 0.09, 0.1, 0.8, 0.81, 0.79, 0.8)
  else
    c(0.1, 0.4, 0.2, 0.8, 1.2, 1.5, 0.9, 1.8)
  seTE <- if (low.heterogeneity) rep(0.2, 8)
  else c(0.1, 0.2, 0.15, 0.25, 0.1, 0.2, 0.15, 0.25)
  metagen(TE, seTE, subgroup = rep(c("A", "B"), each = 4),
            common = FALSE, method.tau = "REML", method.tau.ci = "",
            method.random.ci = methods, adhoc.hakn.ci = adhoc,
            level.ma = level, tau.common = tau.common)
}

expect_subgroup_q_from_stored_se <- function(x) {
  se <- x$seTE.random.w
  if (!is.matrix(se))
    se <- matrix(se, ncol = 1)
  expect_length(x$Q.b.random, ncol(se))
  expect_length(x$pval.Q.b.random, ncol(se))
  dfs <- if (is.list(x$df.Q.b.random)) x$df.Q.b.random else list(x$df.Q.b.random)
  expect_length(dfs, ncol(se))
  for (i in seq_len(ncol(se))) {
    valid <- is.finite(x$TE.random.w) & is.finite(se[, i]) & se[, i] > 0
    expect_true(all(valid))
    weights <- 1 / se[valid, i]^2
    effects <- x$TE.random.w[valid]
    center <- sum(weights * effects) / sum(weights)
    Q <- sum(weights * (effects - center)^2)
    df <- sum(valid) - 1
    expect_equal(unname(x$Q.b.random[i]), unname(Q), tolerance = 1e-10)
    expect_equal(unname(dfs[[i]]), df)
    expect_equal(unname(x$pval.Q.b.random[i]),
                  unname(pchisq(Q, df, lower.tail = FALSE)), tolerance = 1e-10)
  }
}

test_that("subgroup tests agree for methods requested alone, together, or in reverse order", {
  methods <- c("classic", "HK", "KR", "CR0", "CR1")
  for (common.tau in c(FALSE, TRUE)) {
    forward <- subgroup_test_invariance_fixture(methods, tau.common = common.tau)
    reverse <- subgroup_test_invariance_fixture(rev(methods), tau.common = common.tau)
    expect_subgroup_q_from_stored_se(forward)
    expect_subgroup_q_from_stored_se(reverse)
    expect_identical(unname(forward$method.random.ci), methods)
    for (i in seq_along(methods)) {
      single <- subgroup_test_invariance_fixture(methods[i], tau.common = common.tau)
      expect_subgroup_q_from_stored_se(single)
      j <- match(methods[i], reverse$method.random.ci)
      expect_equal(unname(forward$Q.b.random[i]), unname(single$Q.b.random),
                    tolerance = 1e-10)
      expect_equal(unname(reverse$Q.b.random[j]), unname(single$Q.b.random),
                    tolerance = 1e-10)
      expect_equal(unname(forward$pval.Q.b.random[i]), unname(single$pval.Q.b.random),
                    tolerance = 1e-10)
      expect_equal(unname(reverse$pval.Q.b.random[j]), unname(single$pval.Q.b.random),
                    tolerance = 1e-10)
    }
  }
})

test_that("changing only the confidence level does not change ordinary subgroup tests", {
  for (common.tau in c(FALSE, TRUE)) {
    reference <- subgroup_test_invariance_fixture(c("classic", "HK"),
                                                  tau.common = common.tau)
    for (level in c(0.90, 0.95, 0.99)) {
      for (methods in list(c("classic", "HK"), c("HK", "classic"))) {
        ma <- subgroup_test_invariance_fixture(methods, level, common.tau)
        expect_subgroup_q_from_stored_se(ma)
        order <- match(ma$method.random.ci, reference$method.random.ci)
        expect_equal(unname(ma$seTE.random.w),
                      unname(reference$seTE.random.w[, order]), tolerance = 1e-10)
        expect_equal(unname(ma$Q.b.random), unname(reference$Q.b.random[order]),
                      tolerance = 1e-10)
        expect_equal(unname(ma$pval.Q.b.random),
                      unname(reference$pval.Q.b.random[order]), tolerance = 1e-10)
        if (level != 0.95)
          expect_false(isTRUE(all.equal(unname(ma$lower.random.w),
                                         unname(reference$lower.random.w[, order]))))
      }
    }
  }
})

test_that("subgroup tests use the SE for each Hartung-Knapp ad hoc correction", {
  methods <- c("classic", "HK", "HK", "HK", "HK")
  adhoc <- c("", "", "se", "ci", "IQWiG6")
  for (common.tau in c(FALSE, TRUE)) {
    reference <- subgroup_test_invariance_fixture(
      methods, tau.common = common.tau, adhoc = adhoc, low.heterogeneity = TRUE)
    expect_subgroup_q_from_stored_se(reference)
    # The corrections must actually alter the SE in this fixture.
    expect_true(all(reference$seTE.random.w[, 3] > reference$seTE.random.w[, 2]))
    expect_gt(abs(reference$Q.b.random[2] - reference$Q.b.random[3]), 1)
    for (order in list(seq_along(methods), rev(seq_along(methods)))) {
      ma <- subgroup_test_invariance_fixture(
        methods[order], tau.common = common.tau, adhoc = adhoc[order],
        low.heterogeneity = TRUE)
      expect_subgroup_q_from_stored_se(ma)
      expect_equal(unname(ma$Q.b.random), unname(reference$Q.b.random[order]),
                    tolerance = 1e-10)
      for (i in seq_along(order)) {
        single <- subgroup_test_invariance_fixture(
          methods[order[i]], tau.common = common.tau, adhoc = adhoc[order[i]],
          low.heterogeneity = TRUE)
        expect_equal(unname(ma$seTE.random.w[, i]), unname(as.vector(single$seTE.random.w)),
                      tolerance = 1e-10)
        expect_equal(unname(ma$Q.b.random[i]), unname(single$Q.b.random),
                      tolerance = 1e-10)
        expect_equal(unname(ma$pval.Q.b.random[i]), unname(single$pval.Q.b.random),
                      tolerance = 1e-10)
      }
    }
  }
})

test_that("CR2 subgroup tests also agree between single and multiple method requests", {
  skip_if_not_installed("clubSandwich")
  for (common.tau in c(FALSE, TRUE)) {
    single <- subgroup_test_invariance_fixture("CR2", tau.common = common.tau)
    expect_subgroup_q_from_stored_se(single)
    for (methods in list(c("classic", "HK", "CR2"), c("CR2", "HK", "classic"))) {
      ma <- subgroup_test_invariance_fixture(methods, tau.common = common.tau)
      expect_subgroup_q_from_stored_se(ma)
      i <- match("CR2", ma$method.random.ci)
      expect_equal(unname(ma$Q.b.random[i]), unname(single$Q.b.random), tolerance = 1e-10)
      expect_equal(unname(ma$pval.Q.b.random[i]), unname(single$pval.Q.b.random),
                    tolerance = 1e-10)
    }
  }
})
