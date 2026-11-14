test_that("component interval endpoints are paired without inventing a total interval", {
  ma <- metabind_component_fixture()[[2]]
  dat <- metabind(ma, ma)$overall[1, , drop = FALSE]
  # Synthetic endpoints exercise formatting only, not statistical CI coverage.
  dat$lower.tau2 <- list(c(first = 0.01, second = 0.04))
  dat$upper.tau2 <- list(c(first = 0.09, second = 0.16))
  dat$lower.tau <- list(c(first = 0.1, second = 0.2))
  dat$upper.tau <- list(c(first = 0.3, second = 0.4))
  before <- dat
  expect_identical(meta:::format_metabind_variance(
    dat$lower.tau2, upper = dat$upper.tau2, digits = 2),
    "[0.01; 0.09], [0.04; 0.16]")
  flags <- list(Q = FALSE, tau2 = TRUE, tau2.ci = TRUE,
                tau = TRUE, tau.ci = TRUE, I2 = FALSE, H = FALSE)
  formatted <- meta:::format_metabind_heterogeneity(
    dat, flags, digits = list(tau2 = 2, tau = 2, I2 = 2),
    text.tau2 = "Variance", text.tau = "SD")
  intervals <- which(names(formatted) == "95% CI")
  expect_length(intervals, 2L)
  expect_identical(formatted[[intervals[1]]], "[0.01; 0.09], [0.04; 0.16]")
  expect_identical(formatted[[intervals[2]]], "[0.10; 0.30], [0.20; 0.40]")
  expect_identical(dat, before)
})

test_that("raw component differences prevent sharing even when totals agree", {
  ma <- metabind_component_fixture(TRUE)[[2]]
  x <- metabind(ma, ma)
  expect_true(meta:::prepare_metabind_subgroups(x)$shared.heterogeneity)
  x$overall$tau2 <- list(c(between = 0.1, within = 0.3),
                         c(between = 0.3, within = 0.1))
  x$overall$tau <- lapply(x$overall$tau2, sqrt)
  expect_identical(sum(x$overall$tau2[[1]]), sum(x$overall$tau2[[2]]))
  prepared <- meta:::prepare_metabind_subgroups(x)
  expect_true(prepared$shared.overall)
  expect_false(prepared$shared.heterogeneity)

  # Differences hidden by display rounding must also remain distinct.
  x$overall$tau2[[2]] <- x$overall$tau2[[1]] + c(1e-8, -1e-8)
  x$overall$tau <- lapply(x$overall$tau2, sqrt)
  expect_identical(meta:::format_metabind_variance(x$overall$tau2, digits = 4),
                   rep("0.1000, 0.3000", 2))
  expect_false(meta:::prepare_metabind_subgroups(x)$shared.heterogeneity)
})

test_that("list all-NA cleanup retains explicit NULL transformation arguments", {
  dat <- data.frame(id = 1:2)
  dat$tau2 <- list(c(a = NA_real_, b = NA_real_), NA_real_)
  dat$tau <- list(c(a = NA_real_, b = 0.2), NA_real_)
  dat$args.transf <- list(NULL, NULL)
  dat$args.backtransf <- list(NULL, list(base = NULL))
  cleaned <- meta:::drop_allNAs(dat)
  expect_named(cleaned, c("id", "tau", "args.transf", "args.backtransf"))
  expect_identical(cleaned$args.transf, dat$args.transf)
  expect_identical(cleaned$args.backtransf, dat$args.backtransf)
  expect_identical(cleaned$tau, dat$tau)
  expect_identical(meta:::drop_allNAs(dat, keep = "tau2"), dat)

  ma <- metabind_component_fixture()[[1]]
  x <- metabind(ma, ma)
  exported <- as.data.frame(x)
  expect_identical(exported$args.transf, rep(list(NULL), nrow(exported)))
  expect_identical(exported$args.backtransf, rep(list(NULL), nrow(exported)))
})
