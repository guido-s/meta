metabind_heterogeneity_fixture <- function() {
  data.frame(
    id = 1:2, name = c("First", "Second"), func = "metagen", sm = "MD",
    method.tau = "REML", method.tau.ci = "QP", method.I2 = "Q",
    level.hetstat = 0.95, Q = 3.68, df.Q = 4, pval.Q = 0.4514,
    tau2 = 0, lower.tau2 = 0, upper.tau2 = 0.7255, se.tau2 = NA_real_,
    tau = 0, lower.tau = 0, upper.tau = 0.8518,
    H = 1, lower.H = 1, upper.H = 2.19,
    I2 = 0, lower.I2 = 0, upper.I2 = 0.792)
}

format_metabind_test_heterogeneity <- function(dat, print = list(), ...) {
  flags <- list(func = NULL, sm = NULL, method = NULL, level.hetstat = NULL,
                Q = TRUE, pval.Q = TRUE, tau2 = TRUE, tau2.ci = TRUE,
                tau = TRUE, tau.ci = TRUE, I2 = TRUE, I2.ci = TRUE,
                H = TRUE, H.ci = TRUE)
  flags <- utils::modifyList(flags, print, keep.null = TRUE)
  meta:::format_metabind_heterogeneity(
    dat, flags, digits = list(Q = 2, pval.Q = 4, tau2 = 4, tau = 4, I2 = 2), ...)
}

test_that("heterogeneity formatting handles IDs, labels, methods, and repeated CI headings", {
  dat <- metabind_heterogeneity_fixture()
  formatted <- format_metabind_test_heterogeneity(
    dat, ids = 2, print = list(method = TRUE),
    text.tau2 = "Variance", text.tau = "SD", text.I2 = "I-squared")
  expect_identical(formatted$name, "Second")
  expect_true(all(c("Variance", "SD", "I-squared", "H", "Q", "df", "p-value",
                    "method.tau", "method.tau.ci", "method.I2") %in% names(formatted)))
  expect_equal(sum(names(formatted) == "95% CI"), 4)
  expect_false(any(c("func", "sm", "level", "se.tau2", "upper.tau2") %in% names(formatted)))
})

test_that("pooled heterogeneity deduplicates raw fields without collapsing rounded differences", {
  dat <- metabind_heterogeneity_fixture()
  expect_equal(nrow(format_metabind_test_heterogeneity(dat)), 2)
  pooled <- format_metabind_test_heterogeneity(dat, pooled = TRUE)
  expect_equal(nrow(pooled), 1)
  expect_false("name" %in% names(pooled))

  dat$Q[2] <- dat$Q[2] + 1e-8
  pooled <- format_metabind_test_heterogeneity(dat, pooled = TRUE)
  expect_equal(nrow(pooled), 2)
  expect_identical(pooled$Q[1], pooled$Q[2])

  dat <- metabind_heterogeneity_fixture()
  dat$method.tau[2] <- "DL"
  pooled <- format_metabind_test_heterogeneity(
    dat, pooled = TRUE, print = list(method = FALSE))
  expect_equal(nrow(pooled), 2)
  expect_false("method.tau" %in% names(pooled))
})

test_that("heterogeneity formatting respects level and statistic controls", {
  dat <- metabind_heterogeneity_fixture()
  dat$level.hetstat[2] <- 0.90
  formatted <- format_metabind_test_heterogeneity(dat)
  expect_identical(trimws(formatted$level), c("0.9500", "0.9000"))
  expect_equal(sum(names(formatted) == "CI"), 4)

  flags <- list(Q = FALSE, tau2 = FALSE, tau = FALSE, I2 = FALSE, H.ci = FALSE)
  formatted <- format_metabind_test_heterogeneity(dat, print = flags)
  expect_identical(names(formatted), c("name", "H"))

  dat$Q <- dat$df.Q <- dat$pval.Q <- NA_real_
  dat$lower.tau2 <- dat$upper.tau2 <- NA_real_
  formatted <- format_metabind_test_heterogeneity(dat, print = list(method = TRUE))
  expect_false(any(c("Q", "df", "p-value") %in% names(formatted)))
  expect_equal(sum(names(formatted) == "CI"), 3)
})

test_that("heterogeneity degrees of freedom and levels honor missing labels", {
  dat <- metabind_heterogeneity_fixture()
  dat$df.Q <- c(1000, NA_real_)
  dat$level.hetstat <- c(0.95, NA_real_)
  formatted <- format_metabind_test_heterogeneity(
    dat, print = list(level.hetstat = TRUE), lab.NA = "MISSING", big.mark = ",")
  expect_identical(trimws(formatted$df), c("1,000", "MISSING"))
  expect_identical(trimws(formatted$level), c("0.9500", "MISSING"))

  dat$df.Q <- dat$level.hetstat <- NA_real_
  formatted <- format_metabind_test_heterogeneity(
    dat, print = list(level.hetstat = TRUE), lab.NA = "MISSING")
  expect_false(any(c("df", "level") %in% names(formatted)))
})
