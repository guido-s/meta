test_that("metabind heterogeneity controls validate before any value is ignored", {
  grouped <- metabind_render_fixture()
  ma <- metagen_forest_options()
  plain <- metabind(ma, ma)
  invalid <- list("common", "random", "study", NA, logical(), c(TRUE, FALSE), list(TRUE))
  for (x in list(grouped, plain)) {
    for (value in c(list(NULL), invalid))
      expect_error(forest(x, hetstat = value), "Argument 'hetstat'", fixed = TRUE)
    for (value in c(list(NULL), invalid))
      expect_error(forest(x, overall = value), "Argument 'overall'", fixed = TRUE)
    for (value in invalid)
      expect_error(forest(x, overall.hetstat = value),
                   "Argument 'overall.hetstat'", fixed = TRUE)
  }
})

test_that("shared overall heterogeneity requires both logical controls", {
  x <- metabind_render_fixture()
  settings <- list(list(), list(overall.hetstat = NULL),
                   list(overall.hetstat = TRUE), list(overall.hetstat = FALSE),
                   list(hetstat = FALSE), list(hetstat = FALSE, overall.hetstat = TRUE))
  visible <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  for (i in seq_along(settings)) {
    expect_silent(drawn <- do.call(
      metabind_render_pdf, c(list(x = x, details = FALSE), settings[[i]])))
    expect_identical(drawn$args$overall.hetstat, visible[i])
    expect_false(drawn$args$subgroup.hetstat)
    labels <- vapply(drawn$labels[!is.na(drawn$columns$col.studlab$rows)], function(label)
      paste(deparse(label), collapse = " "), character(1))
    expect_identical(any(grepl("Heterogeneity", labels, fixed = TRUE)), visible[i])
  }
})

test_that("overall heterogeneity controls do not disable component heterogeneity", {
  x <- metabind_render_fixture(method.tau = "REML")
  for (overall.hetstat in list(NULL, TRUE, FALSE)) {
    expect_silent(drawn <- metabind_render_pdf(
      x, overall.hetstat = overall.hetstat, details = FALSE))
    expect_true(drawn$args$hetstat)
    expect_false(drawn$args$overall.hetstat)
    expect_true(drawn$args$subgroup.hetstat)
    labels <- vapply(drawn$labels[!is.na(drawn$columns$col.studlab$rows)], function(label)
      paste(deparse(label), collapse = " "), character(1))
    expect_true(any(grepl("Heterogeneity", labels, fixed = TRUE)))
  }
  expect_silent(disabled <- metabind_render_pdf(
    x, hetstat = FALSE, overall.hetstat = TRUE, details = FALSE))
  expect_false(disabled$args$hetstat)
  expect_false(disabled$args$overall.hetstat)
  expect_false(disabled$args$subgroup.hetstat)
})

test_that("matching effects do not share differing heterogeneity even when enabled", {
  x <- metabind_render_fixture(method.tau = "REML")
  drawn <- metabind_render_pdf(
    x, random = FALSE, prediction = FALSE, hetstat = TRUE,
    overall.hetstat = TRUE, details = FALSE)
  expect_true(drawn$args$overall)
  expect_false(drawn$args$hetstat)
  expect_false(drawn$args$overall.hetstat)
  expect_false(drawn$args$subgroup.hetstat)
})

test_that("plain metabind plots ignore row controls but retain requested columns", {
  ma <- metagen_forest_options()
  x <- metabind(ma, ma)
  expect_silent(default <- metabind_render_pdf(x, details = FALSE))
  expect_false(default$args$hetstat)
  expect_false(default$args$overall.hetstat)
  expect_silent(automatic <- metabind_render_pdf(
    x, overall.hetstat = NULL, details = FALSE))
  expect_false(automatic$args$overall.hetstat)
  for (value in c(TRUE, FALSE)) {
    expect_warning(drawn <- metabind_render_pdf(
      x, overall.hetstat = value, details = FALSE),
      "Input to argument 'overall.hetstat' ignored", fixed = TRUE)
    expect_false(drawn$args$overall.hetstat)
  }
  expect_warning(drawn <- metabind_render_pdf(x, hetstat = TRUE, details = FALSE),
                 "Input to argument 'hetstat' ignored", fixed = TRUE)
  expect_false(drawn$args$hetstat)
  expect_warning(drawn <- metabind_render_pdf(x, overall = TRUE, details = FALSE),
                 "Input to argument 'overall' ignored", fixed = TRUE)
  expect_false(drawn$args$overall)
  expect_silent(disabled <- metabind_render_pdf(x, overall = FALSE, details = FALSE))
  expect_false(disabled$args$overall)
  columns <- metabind_render_pdf(
    x, hetstat = FALSE, details = FALSE, common = FALSE, prediction = FALSE,
    print = list(I2 = TRUE), leftcols = c("studlab", "I2"), lab.NA = "MISSING")
  expect_true("I2" %in% columns$result$leftcols)
  expect_true(all(columns$result$data$I2 != "MISSING"))
})

test_that("native forest retains its character-valued heterogeneity placement", {
  ma <- metagen_forest_options()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(forest(ma, hetstat = "random", details = FALSE))
})
