test_that("calcwidth lists retain pooled inheritance and legacy controls", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3), method.tau = "DL")
  observed <- logical()
  calculate <- meta:::calc_pooled_prefix_extra
  testthat::local_mocked_bindings(
    calc_pooled_prefix_extra = function(lines, pooled.cols, do.calc, ...) {
      observed <<- c(observed, do.calc)
      calculate(lines, pooled.cols, do.calc, ...)
    }, .package = "meta")
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())

  settings <- list(list(pooled = FALSE),
                    list(pooled = FALSE, random = TRUE),
                    list(pooled = TRUE, common = FALSE))
  expected <- list(c(FALSE, FALSE), c(FALSE, TRUE), c(FALSE, TRUE))
  for (legacy in c(FALSE, TRUE)) {
    for (i in seq_along(settings)) {
      if (legacy)
        options <- setNames(settings[[i]], paste0("calcwidth.", names(settings[[i]])))
      else
        options <- list(calcwidth = settings[[i]])
      observed <- logical()
      do.call(meta::forest, c(list(x = ma), options))
      expect_identical(observed, expected[[i]])
    }
  }
  observed <- logical()
  forest(ma, common = FALSE, calcwidth = list(pooled = TRUE))
  expect_identical(observed, c(FALSE, TRUE))
  expect_error(forest(ma, calcwidth = list(common = "invalid")), "calcwidth[.]common")
})

forest_calcwidth_measure <- function(x, ...) {
  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off())
  result <- meta::forest(x, details = FALSE, ...)
  plot <- grid::grid.grab()
  layout <- plot$childrenvp[[1]]$parent$layout
  labels <- as.list(plot$children)
  heterogeneity <- vapply(labels, function(grob)
    is.call(grob$label) && is.character(grob$label[[2]]) &&
      startsWith(grob$label[[2]], "Heterogeneity"), logical(1))
  list(
    left = grid::convertWidth(
      sum(layout$widths[seq_len(2 * length(result$leftcols) - 1)]),
      "mm", valueOnly = TRUE),
    heterogeneity = vapply(labels[heterogeneity], function(grob)
      grid::convertWidth(grid::grobWidth(grob), "mm", valueOnly = TRUE), numeric(1)))
}

test_that("RevMan5 heterogeneity width includes I2 confidence intervals", {
  ma <- metaprop(c(13, 47, 15, 33, 35), c(129, 400, 103, 177, 38))
  for (totals in c(TRUE, FALSE)) {
    for (legacy in c(FALSE, TRUE)) {
      options <- if (legacy)
        list(calcwidth.hetstat = TRUE, calcwidth.tests = FALSE)
      else
        list(calcwidth = list(hetstat = TRUE, tests = FALSE))
      base <- list(x = ma, layout = "RevMan5", pooled.totals = totals)
      without <- do.call(forest_calcwidth_measure,
                         c(base, list(print.I2.ci = FALSE), options))
      with <- do.call(forest_calcwidth_measure,
                      c(base, list(print.I2.ci = TRUE), options))
      expect_length(with$heterogeneity, 1)
      expect_gt(with$heterogeneity, without$heterogeneity)
      expect_gt(with$left, without$left)
      expect_gte(with$left + 1e-7, with$heterogeneity)
    }

    without <- forest_calcwidth_measure(
      ma, layout = "RevMan5", pooled.totals = totals, print.I2.ci = FALSE)
    with <- forest_calcwidth_measure(
      ma, layout = "RevMan5", pooled.totals = totals, print.I2.ci = TRUE)
    expect_gt(with$left, without$left)

    without <- forest_calcwidth_measure(
      ma, layout = "RevMan5", pooled.totals = totals, print.I2.ci = FALSE,
      calcwidth = list(hetstat = FALSE, tests = FALSE))
    with <- forest_calcwidth_measure(
      ma, layout = "RevMan5", pooled.totals = totals, print.I2.ci = TRUE,
      calcwidth = list(hetstat = FALSE, tests = FALSE))
    expect_equal(with$left, without$left, tolerance = 1e-7)
  }
})

test_that("heterogeneity and test width controls select their own labels in each layout", {
  ma <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), common = FALSE,
                  method.tau = "DL", null.effect = 0)
  for (layout in c("RevMan5", "meta", "meta8", "BMJ", "JAMA")) {
    for (measure in c("hetstat", "tests")) {
      options <- list(hetstat = FALSE, tests = FALSE)
      options[[measure]] <- TRUE
      short <- forest_calcwidth_measure(
        ma, layout = layout, calcwidth = options, test.overall.random = TRUE,
        label.test.overall.random = "TEST", hetlab = "HET")
      long.test <- forest_calcwidth_measure(
        ma, layout = layout, calcwidth = options, test.overall.random = TRUE,
        label.test.overall.random = strrep("LONG_TEST_", 15), hetlab = "HET")
      long.het <- forest_calcwidth_measure(
        ma, layout = layout, calcwidth = options, test.overall.random = TRUE,
        label.test.overall.random = "TEST", hetlab = strrep("LONG_HET_", 15))
      if (measure == "hetstat") {
        expect_equal(long.test$left, short$left, tolerance = 1e-7)
        expect_gt(long.het$left, short$left)
      }
      else {
        expect_gt(long.test$left, short$left)
        expect_equal(long.het$left, short$left, tolerance = 1e-7)
      }
    }
  }
})
