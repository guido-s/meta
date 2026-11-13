three_level_forest_fixture <- function(methods = "HK", rho = 0.3,
                                        three.level = TRUE) {
  dat <- data.frame(
    TE = c(0.1, 0.4, 0.2, 0.8, 0.3, 0.6, 1.2, 1.5, 0.9, 1.8, 1.1, 1.6),
    seTE = rep(c(0.1, 0.2, 0.15), 4), cluster = rep(1:6, each = 2),
    group = rep(c("B", "A"), each = 6))
  metagen(TE, seTE, data = dat, sm = "MD",
          cluster = if (three.level) cluster else NULL, rho = rho,
          subgroup = group, tau.common = TRUE, method.tau.ci = "",
          method.random.ci = methods)
}

three_level_forest_capture <- function(x, fields, ...) {
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      args <- list(...)
      observed <<- list(values = mget(fields, envir = parent.frame()), args = args)
      do.call(draw, args)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 20, height = 16)
  on.exit(grDevices::dev.off(), add = TRUE)
  forest(x, ...)
  observed
}

test_that("three-level forest annotations use total overall and residual variance", {
  x <- three_level_forest_fixture()
  before <- x
  expect_length(x$tau2, 2L)
  expect_length(x$tau2.resid, 2L)
  fields <- c("tau2", "tau", "tau2.resid", "tau.resid",
              "lower.tau2", "upper.tau2", "lower.tau", "upper.tau",
              "lower.tau2.resid", "upper.tau2.resid",
              "lower.tau.resid", "upper.tau.resid", "hetstat.resid")
  # Component-wise endpoints must not become an aggregate confidence interval.
  for (field in c("tau2", "tau2.resid", "tau", "tau.resid")) {
    x[[paste0("lower.", field)]] <- c(0.01, 0.02)
    x[[paste0("upper.", field)]] <- c(0.8, 0.9)
  }
  input <- x
  for (layout in c("meta", "RevMan5", "JAMA")) for (ci in c("tau2", "tau")) {
    expect_silent(drawn <- three_level_forest_capture(
      x, fields, layout = layout, details = TRUE, resid.hetstat = TRUE,
      print = list(tau2 = TRUE, tau = TRUE,
                   tau2.ci = ci == "tau2", tau.ci = ci == "tau")))
    values <- drawn$values
    expect_equal(values$tau2, sum(x$tau2))
    expect_equal(values$tau, sqrt(sum(x$tau^2)))
    expect_equal(values$tau2.resid, sum(x$tau2.resid))
    expect_equal(values$tau.resid, sqrt(sum(x$tau.resid^2)))
    for (field in grep("^lower\\.|^upper\\.", fields, value = TRUE))
      expect_true(all(is.na(values[[field]])), info = field)
    expect_false(identical(values$hetstat.resid, ""))
  }
  expect_identical(x, input)
  expect_silent(hidden <- three_level_forest_capture(
    before, "hetstat.resid", resid.hetstat = FALSE))
  expect_identical(hidden$values$hetstat.resid, "")
})

test_that("scalar tau confidence intervals use their own upper-bound signs", {
  x <- three_level_forest_fixture(three.level = FALSE)
  x$lower.tau2 <- 0.01
  x$upper.tau2 <- 0.81
  x$lower.tau <- 0.1
  x$upper.tau <- 0.9
  x$lower.tau2.resid <- 0.04
  x$upper.tau2.resid <- 0.64
  x$lower.tau.resid <- 0.2
  x$upper.tau.resid <- 0.8
  x$sign.lower.tau <- x$sign.lower.tau.resid <- ""
  x$sign.upper.tau <- ">"
  x$sign.upper.tau.resid <- "<"
  before <- x
  for (ci in c("tau2", "tau")) {
    fields <- paste0("hetstat.", ci, c("", ".resid"))
    drawn <- three_level_forest_capture(
      x, fields, resid.hetstat = TRUE, print = list(tau2 = TRUE, tau = TRUE,
        tau2.ci = ci == "tau2", tau.ci = ci == "tau"))
    expect_match(drawn$values[[fields[1]]],
                 if (ci == "tau2") "; >0.8100]" else "; >0.9000]", fixed = TRUE)
    expect_match(drawn$values[[fields[2]]],
                 if (ci == "tau2") "; <0.6400]" else "; <0.8000]", fixed = TRUE)
  }
  expect_identical(x, before)
})

test_that("merged forest subgroup values follow source and CI-method slots", {
  for (three.level in c(FALSE, TRUE)) {
    single <- three_level_forest_fixture("HK", three.level = three.level)
    multiple <- three_level_forest_fixture(c("HK", "classic"), rho = 0.7,
                                            three.level = three.level)
    if (!three.level)
      multiple <- update(multiple, method.tau = "DL")
    reversed <- three_level_forest_fixture(c("classic", "HK"),
                                            three.level = three.level)
    expect_gt(abs(single$TE.random.w[1] - multiple$TE.random.w[1]), 1e-5)
    for (sources in list(list(single, multiple), list(multiple, single),
                         list(multiple, reversed))) {
      x <- metamerge(sources[[1]], sources[[2]], label1 = "Left", label2 = "Right")
      before <- x
      fields <- unlist(lapply(c(if (x$common) "common", "random"), function(model)
        paste0(c("TE", "lower", "upper", "statistic", "pval"), ".", model, ".w")))
      for (sorted in c(FALSE, TRUE)) {
        expect_silent(drawn <- three_level_forest_capture(
          x, c(fields, "subgroup.levels"), sort.subgroup = sorted))
        levels <- drawn$values$subgroup.levels
        for (field in fields) {
          model <- if (grepl("common", field)) "common" else "random"
          n.method <- vapply(sources, function(source)
            length(source[[paste0("lower.", model)]]), integer(1))
          blocks <- lapply(seq_along(sources), function(i) {
            value <- sources[[i]][[field]]
            if (is.matrix(value))
              value[match(levels, rownames(value)), , drop = FALSE]
            else
              matrix(value[match(levels, names(value))], nrow = length(levels),
                     ncol = n.method[i])
          })
          expected <- as.vector(t(do.call(cbind, blocks)))
          expect_equal(unname(drawn$values[[field]]), expected, info = field)
        }
        # The displayed effect and CI columns retain matching row arrays.
        expect_length(drawn$args$cols$col.effect$labels,
                      length(drawn$args$cols$col.ci$labels))
      }
      expect_identical(x, before)
    }
  }
})

test_that("merged forest placeholders handle hidden models and nested merges", {
  for (three.level in c(FALSE, TRUE)) {
    first <- three_level_forest_fixture(three.level = three.level)
    second <- three_level_forest_fixture(c("classic", "HK"), rho = 0.7,
                                          three.level = three.level)
    x <- metamerge(metamerge(first, second), first)
    before <- x
    for (selection in list(list(common = FALSE), list(random = FALSE),
                           list(common = FALSE, random = FALSE),
                           list(overall = FALSE), list(subgroup = c(TRUE, FALSE)))) {
      expect_silent(drawn <- do.call(
        three_level_forest_capture,
        c(list(x = x, fields = c("TEs", "lowTEs", "uppTEs")), selection)))
      expect_length(drawn$values$TEs, length(drawn$values$lowTEs))
      expect_length(drawn$values$TEs, length(drawn$values$uppTEs))
    }
    expect_identical(x, before)
  }
})
