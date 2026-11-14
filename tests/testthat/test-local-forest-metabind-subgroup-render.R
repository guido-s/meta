test_that("shared summaries describe methods rather than showing redundant columns", {
  x <- metabind_render_fixture()
  x <- metabind(update(x$meta[[1]], method = "MH", method.random.ci = "HK"),
                 update(x$meta[[2]], method = "MH", method.random.ci = "HK"))
  drawn <- metabind_render_pdf(x)
  expect_true(drawn$args$details)
  details <- paste(drawn$args$x$.text.details.methods, collapse = " ")
  expect_match(details, "Mantel-Haenszel")
  expect_match(details, "Inverse variance")
  expect_match(details, "Hartung")
  expect_false(any(c("method", "method.interval", "statistic", "pval", "df") %in%
                     c(drawn$result$leftcols, drawn$result$rightcols)))

  disabled <- metabind_render_pdf(x, details = FALSE)
  expect_false(disabled$args$details)
  expect_null(disabled$args$x$.text.details.methods)
})

test_that("shared summaries retain interval methods varying within a result type", {
  x <- metabind_render_fixture(method.random.ci = c("classic", "HK"))
  drawn <- metabind_render_pdf(x)
  expect_true(drawn$args$details)
  expect_true("method.interval" %in% drawn$result$leftcols)
  expect_false("method" %in% drawn$result$leftcols)
  expect_false(any(c("statistic", "pval", "df") %in% drawn$result$leftcols))
})

test_that("nonshared summaries show effect tests in rows rather than automatic columns", {
  x <- metabind_render_fixture(method.tau = "REML")
  drawn <- metabind_render_pdf(x)
  expect_true(drawn$args$details)
  expect_false(any(c("statistic", "pval") %in%
                     c(drawn$result$leftcols, drawn$result$rightcols)))
  expect_true("df" %in% drawn$result$leftcols)
  labels <- Filter(function(z) is.call(z) &&
                     grepl("Test for overall effect", as.character(z[[2]]), fixed = TRUE),
                   drawn$labels)
  expect_length(labels, 4)
})

test_that("explicit shared-summary display options override automatic selection", {
  x <- metabind_render_fixture()
  drawn <- metabind_render_pdf(
    x, print = list(method = TRUE, method.interval = TRUE, stat = TRUE, pval = TRUE),
    rightcols = "df")
  expect_true(all(c("method", "method.interval", "statistic", "pval") %in%
                    drawn$result$leftcols))
  expect_equal(drawn$result$rightcols, "df")
})

test_that("subgroup test columns respect model selection and explicit columns", {
  x <- metabind_render_fixture()
  for (common in c(FALSE, TRUE)) {
    drawn <- metabind_render_pdf(x, common = common, random = !common,
                                 prediction = FALSE)
    expect_true(all(c("Q.b", "df.Q.b", "pval.Q.b") %in% drawn$result$leftcols))
    first <- as.character(drawn$result$data$subgroup.levels) %in% c("A", "C")
    expect_true(all(nzchar(drawn$result$data$Q.b[first])))
    expect_true(all(drawn$result$data$Q.b[!first] == ""))
  }
  pi <- metabind_render_pdf(x, common = FALSE, random = FALSE)
  expect_false(any(c("Q.b", "df.Q.b", "pval.Q.b") %in% pi$result$leftcols))
  hidden <- metabind_render_pdf(x, print = list(pval.Q.b = FALSE))
  expect_true(all(c("Q.b", "df.Q.b") %in% hidden$result$leftcols))
  expect_false("pval.Q.b" %in% hidden$result$leftcols)
  selected <- metabind_render_pdf(x, leftcols = c("studlab", "effect", "ci"))
  expect_equal(selected$result$leftcols, c("studlab", "effect", "ci"))
})

test_that("subgroup labels identify result types without a separate type column", {
  x <- metabind_render_fixture()
  for (grouped in c(FALSE, TRUE)) {
    if (grouped)
      x$samedata <- FALSE
    drawn <- metabind_render_pdf(x)
    dat <- drawn$result$data
    expect_false("type" %in% drawn$result$leftcols)
    expect_false("type" %in% drawn$result$rightcols)
    for (i in c("common", "random", "prediction")) {
      model <- switch(i, common = "Common effect model",
                      random = "Random effects model", prediction = "Prediction interval")
      rows <- dat$type == model
      expect_equal(drawn$args$x$studlab[rows],
                   paste0(dat$subgroup.levels[rows], " (", i, ")"))
    }
  }
  explicit <- metabind_render_pdf(x, studlab = subgroup.levels,
                                  print = list(type = TRUE))
  expect_equal(explicit$args$x$studlab,
               as.character(explicit$result$data$subgroup.levels))
  expect_true("type" %in% explicit$result$leftcols)
  for (i in c("common", "random", "prediction")) {
    single <- metabind_render_pdf(x, common = i == "common",
                                  random = i == "random", prediction = i == "prediction")
    expect_equal(single$args$x$studlab,
                 as.character(single$result$data$subgroup.levels))
    expect_false("type" %in% single$result$leftcols)
  }
})

test_that("effect test rows follow model selections and omit unavailable components", {
  x <- metabind_render_fixture()
  x$samedata <- FALSE
  x$common <- subset(x$common, id == 1)
  x$subgroup$common <- subset(x$subgroup$common, id == 1)
  drawn <- metabind_render_pdf(x)
  expect_equal(drawn$args$test$effect.subgroup.common, c(TRUE, FALSE))
  expect_equal(drawn$args$test$effect.subgroup.random, c(TRUE, TRUE))
  for (i in c("common", "random", "prediction")) {
    drawn <- metabind_render_pdf(x, common = i == "common",
                                 random = i == "random", prediction = i == "prediction")
    labels <- Filter(function(z) is.call(z) &&
                       grepl("Test for overall effect", as.character(z[[2]]), fixed = TRUE),
                     drawn$labels)
    expect_length(labels, switch(i, common = 1, random = 2, prediction = 0))
  }
})

test_that("subgroup heterogeneity uses the requested forest defaults and overrides", {
  x <- metabind_render_fixture()
  drawn <- metabind_render_pdf(x)
  expect_equal(drawn$args$print$I2, meta:::replaceNULL(gs("print.I2"), TRUE))
  expect_true(drawn$args$print$pval.Q)
  expect_false(drawn$args$print$Q)

  old.I2 <- gs("print.I2")
  withr::defer(meta:::setOption("print.I2", old.I2))
  meta:::setOption("print.I2", FALSE)
  expect_false(metabind_render_pdf(x)$args$print$I2)
  override <- metabind_render_pdf(x, print = list(I2 = TRUE, pval.Q = FALSE))
  expect_true(override$args$print$I2)
  expect_false(override$args$print$pval.Q)
  nullable <- metabind_render_pdf(x, print = list(I2 = NULL, pval.Q = NULL))
  expect_true(nullable$args$print$I2)
  expect_true(nullable$args$print$pval.Q)
  meta:::setOption("print.I2", NULL)
  expect_true(metabind_render_pdf(x)$args$print$I2)
})

test_that("overall-effect controls also control component tests and null details", {
  controls <- list(
    list(test.overall = FALSE),
    list(test.overall.common = FALSE),
    list(test.overall.random = FALSE),
    list(test.overall = FALSE, test.overall.common = TRUE))
  for (grouped in c(FALSE, TRUE)) {
    x <- metabind_render_fixture()
    if (grouped)
      x$samedata <- FALSE
    for (i in seq_along(controls)) {
      drawn <- do.call(metabind_render_pdf, c(list(x = x), controls[[i]]))
      labels <- Filter(function(z) is.call(z) &&
        grepl("Test for overall effect", as.character(z[[2]]), fixed = TRUE), drawn$labels)
      expect_length(labels, if (i == 1) 0 else if (grouped) 2 else 1)
      if (i != 1) {
        family <- if (i == 2) "random" else "common"
        expect_true(all(vapply(labels, function(z)
          grepl(family, z[[2]], fixed = TRUE), logical(1))))
      }
      text <- paste(drawn$args$x$.text.details.methods, collapse = " ")
      expect_identical(grepl("Null hypothesis", text, fixed = TRUE), i != 1)
      expect_true(drawn$args$details)
    }
    disabled <- metabind_render_pdf(x, details = FALSE)
    expect_false(disabled$args$details)
    expect_null(disabled$args$x$.text.details.methods)
  }
  x <- metabind_render_fixture()
  drawn <- metabind_render_pdf(x, overall = FALSE)
  expect_false(grepl("Null hypothesis",
                     paste(drawn$args$x$.text.details.methods, collapse = " "), fixed = TRUE))
})
