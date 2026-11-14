forest_label_fixture <- function() {
  metabin(c(10, 20, 15, 25), rep(100, 4),
           c(8, 12, 10, 15), rep(100, 4),
           subgroup = c("A", "A", "B", "B"),
           method = "Inverse", method.tau = "DL")
}

forest_label_plot <- function(x, ...) {
  observed <- axis.input <- NULL
  draw <- meta:::forest_meta_internal
  axis <- meta:::axis_labels
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      observed <<- list(...)
      do.call(draw, observed)
    },
    axis_labels = function(xlim, at, label, ...) {
      axis.input <<- label
      axis(xlim, at, label, ...)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 16, height = 12)
  on.exit(grDevices::dev.off())
  result <- forest(x, ...)
  list(result = result, args = observed, axis = axis.input)
}

test_that("label lists and legacy label fields preserve explicit headings and tests", {
  x <- forest_label_fixture()
  settings <- list(e = "Treatment", c = "Comparator", n = "Participants",
                    e.attach = "n.e", c.attach = "n.c",
                    left = "Favours treatment", right = "Favours comparator",
                    test.overall.common = "COMMON_TEST",
                    test.overall.random = "RANDOM_TEST")
  for (legacy in c(FALSE, TRUE)) {
    if (legacy)
      options <- setNames(settings, paste0("label.", names(settings)))
    else
      options <- list(label = settings)
    drawn <- do.call(forest_label_plot,
                     c(list(x = x, leftcols = c("studlab", "n.e", "n.c", "effect", "ci"),
                            test = list(overall = TRUE)), options))
    expect_true(drawn$axis)
    expect_equal(drawn$result$leftlabs[2:3], rep("Participants", 2))
    expect_identical(drawn$args$col.label.e$labels[[1]]$label, "Treatment")
    expect_identical(drawn$args$col.label.c$labels[[1]]$label, "Comparator")
    labels <- lapply(drawn$args$cols$col.studlab$labels, `[[`, "label")
    tags <- vapply(Filter(is.call, labels), function(z) as.character(z[[2]]), "")
    expect_true(all(c("COMMON_TEST", "RANDOM_TEST") %in% tags))
  }
})

test_that("logical and custom axis labels remain independent of label settings", {
  x <- forest_label_fixture()
  expect_true(forest_label_plot(x)$axis)
  expect_false(forest_label_plot(x, label = FALSE)$axis)
  expect_warning(old <- forest_label_plot(x, labels = FALSE), "instead of 'labels'")
  expect_false(old$axis)
  expect_warning(modern <- forest_label_plot(x, label = list(n = "Participants"),
                                              labels = FALSE), "ignored")
  expect_true(modern$axis)
  ticks <- c("Low", "Middle", "High")
  modern <- forest_label_plot(x, label = ticks, at = c(-1, 0, 1),
                              xlim = c(-1, 1), backtransf = FALSE)
  expect_identical(modern$axis, ticks)
  expect_warning(old <- forest_label_plot(x, labels = ticks, at = c(-1, 0, 1),
                                           xlim = c(-1, 1), backtransf = FALSE),
                 "instead of 'labels'")
  expect_identical(old$axis, ticks)
})

test_that("deprecated treatment labels are resolved before layout processing", {
  x <- forest_label_fixture()
  drawn <- forest_label_plot(x, lab.e = "Legacy treatment", lab.c = "Legacy comparator",
                              lab.e.attach.to.col = "n.e", lab.c.attach.to.col = "n.c",
                              leftcols = c("studlab", "n.e", "n.c", "effect", "ci"),
                              warn.deprecated = FALSE)
  expect_identical(drawn$args$col.label.e$labels[[1]]$label, "Legacy treatment")
  expect_identical(drawn$args$col.label.c$labels[[1]]$label, "Legacy comparator")
  expect_warning(drawn <- forest_label_plot(x, label = list(e = "Modern treatment"),
                                            lab.e = "Legacy treatment"), "ignored")
  expect_identical(drawn$args$col.label.e$labels[[1]]$label, "Modern treatment")
  expect_error(forest_label_plot(x, label = list(e = "Modern"), label.e = "Legacy"),
               "cannot both be specified")
})

test_that("generated test labels use the resolved model selection", {
  x <- forest_label_fixture()
  for (both in c(FALSE, TRUE)) {
    drawn <- forest_label_plot(x, random = both, test = list(overall = TRUE))
    labels <- lapply(drawn$args$cols$col.studlab$labels, `[[`, "label")
    labels <- Filter(function(z) is.call(z) &&
                       grepl("Test for overall effect", as.character(z[[2]]), fixed = TRUE),
                     labels)
    expect_length(labels, if (both) 2 else 1)
    expect_identical(grepl("(common effect)", labels[[1]][[2]], fixed = TRUE), both)
  }
})

test_that("short heading vectors label additional columns on either side", {
  data(Olkin1995, package = "meta")
  x <- metabin(ev.exp, n.exp, ev.cont, n.cont,
                data = Olkin1995, subset = c(41, 47, 51, 59),
                sm = "RR", method = "I", studlab = paste(author, year))
  newcols <- c("author", "year")
  #
  for (side in c("left", "right")) {
    options <- if (side == "left") list() else
      list(leftcols = "studlab", rightcols = c("effect", "ci"))
    baseline <- do.call(forest_label_plot, c(list(x = x), options))
    cols <- paste0(side, "cols")
    labs <- paste0(side, "labs")
    options[[cols]] <- c(baseline$result[[cols]], newcols)
    #
    for (missing.label in c(FALSE, TRUE)) {
      options[[labs]] <- c(if (missing.label) NA_character_ else "Author",
                           "Year of publication")
      expected <- c(if (missing.label) "author" else "Author",
                     "Year of publication")
      drawn <- do.call(forest_label_plot, c(list(x = x), options))
      #
      expect_identical(drawn$args[[paste0(labs, ".new")]], expected)
      expect_identical(tail(drawn$result[[labs]], 2), expected)
      expect_identical(head(drawn$result[[labs]], -2), baseline$result[[labs]])
      #
      headings <- vapply(drawn$args$cols[paste0("col.", newcols)],
                          function(z) z$labels[[1]]$label, character(1),
                          USE.NAMES = FALSE)
      expect_identical(headings, expected)
    }
  }
})
