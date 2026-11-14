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
