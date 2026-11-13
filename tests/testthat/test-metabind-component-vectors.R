metabind_component_fields <- c(
  "tau2", "lower.tau2", "upper.tau2", "se.tau2",
  "tau", "lower.tau", "upper.tau")

metabind_component_fixture <- function(subgroup = FALSE) {
  dat <- data.frame(
    TE = c(-0.4, 0.1, 0.2, 0.7, -0.2, 0.4, 0.8, 1.1, 0.3, 0.9, 1.2, 0.6),
    seTE = rep(c(0.15, 0.2), 6), cluster = rep(1:6, each = 2),
    group = rep(c("A", "B"), each = 6))
  args <- list(TE = dat$TE, seTE = dat$seTE, data = dat, sm = "MD",
               common = FALSE, random = TRUE, tau.common = FALSE,
               method.random.ci = "HK", prediction = TRUE, method.predict = "V")
  if (subgroup)
    args$subgroup <- dat$group
  two <- do.call(metagen, args)
  args$cluster <- dat$cluster
  three <- do.call(metagen, args)
  sources <- list(two, three)
  # Name even scalar fields so extraction cannot silently discard names.
  for (i in seq_along(sources))
    for (field in metabind_component_fields) {
      value <- as.numeric(sources[[i]][[field]])
      sources[[i]][[field]] <- setNames(value, paste0("component", seq_along(value)))
    }
  sources
}

test_that("variance list columns retain complete source vectors in every result row", {
  for (subgroup in c(FALSE, TRUE)) {
    sources <- metabind_component_fixture(subgroup)
    x <- metabind(sources[[1]], sources[[2]], keepmeta = TRUE)
    before <- x
    expect_length(sources[[1]]$tau2, 1L)
    expect_length(sources[[2]]$tau2, 2L)
    tables <- list(x$overall, as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE))
    if (subgroup)
      tables <- c(tables, list(as.data.frame(x, subgroup = TRUE, drop_all_NA = FALSE)))
    for (dat in tables)
      for (field in metabind_component_fields) {
        expect_type(dat[[field]], "list")
        for (row in seq_len(nrow(dat)))
          expect_identical(dat[[field]][[row]], sources[[dat$id[row]]][[field]],
                           info = field)
      }
    if (subgroup) {
      prepared <- meta:::prepare_metabind_subgroups(x)
      expect_false(prepared$shared.overall)
      for (dat in list(prepared$dat, prepared$subdat))
        for (field in metabind_component_fields)
          for (row in seq_len(nrow(dat)))
            expect_identical(dat[[field]][[row]], sources[[dat$id[row]]][[field]])
    }
    expect_silent(output <- capture.output(print(
      x, print = list(tau2 = TRUE, tau = TRUE), details.methods = FALSE)))
    expected <- meta:::format_metabind_variance(list(sources[[2]]$tau2), digits = 4)
    if (!subgroup)
      expect_true(any(grepl(expected, output, fixed = TRUE)))
    expect_identical(x, before)
    expect_identical(x$meta, sources)
  }
})

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

test_that("native forest totals are local and have no aggregate component CI", {
  sources <- metabind_component_fixture(TRUE)
  x <- metabind(sources[[1]], sources[[2]], keepmeta = TRUE)
  before <- x
  observed <- NULL
  draw <- meta:::forest.meta
  testthat::local_mocked_bindings(forest.meta = function(...) {
    observed <<- list(...)
    do.call(draw, observed)
  }, .package = "meta")
  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(forest(x, layout = "meta", details = FALSE,
                      print = list(tau2 = TRUE, tau2.ci = TRUE,
                                   tau = TRUE, tau.ci = TRUE)))
  native <- observed$x
  expect_equal(unname(native$tau2.w), vapply(sources, function(ma) sum(ma$tau2), 0))
  expect_equal(unname(native$tau.w),
               vapply(sources, function(ma) sqrt(sum(ma$tau^2)), 0))
  for (field in c("lower.tau2", "upper.tau2", "se.tau2", "lower.tau", "upper.tau")) {
    expect_identical(unname(native[[paste0(field, ".w")]][2]), NA_real_, info = field)
    expect_equal(unname(native[[paste0(field, ".w")]][1]),
                 unname(sources[[1]][[field]]), info = field)
  }
  expect_identical(x, before)
  expect_identical(x$meta, sources)

  sources <- metabind_component_fixture(FALSE)
  x <- metabind(sources[[1]], sources[[2]], keepmeta = TRUE)
  before <- x
  expect_silent(forest(x, layout = "meta", details = FALSE,
                      print = list(tau2 = TRUE, tau2.ci = TRUE,
                                   tau = TRUE, tau.ci = TRUE)))
  expect_identical(x, before)
  expect_identical(x$meta, sources)
})
