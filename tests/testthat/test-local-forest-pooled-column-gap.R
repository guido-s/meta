pooled_gap_fixture <- function(method.random.ci) {
  dat <- data.frame(
    TE = c(-0.4, 0.8, 0.5, -0.3, 0.9, 0.2, -0.6, 0.7),
    seTE = c(0.1, 0.3, 0.2, 0.2, 0.4, 0.15, 0.25, 0.3),
    first = rep(c("A", "B"), each = 4),
    second = rep(c("C", "D"), 4))
  ma <- metagen(TE, seTE, data = dat, subgroup = first,
                method.tau = "DL", method.random.ci = method.random.ci,
                common = TRUE, random = TRUE, prediction = FALSE)
  metabind(ma, update(ma, subgroup = second, method.tau = "REML"),
           name = c("First grouping", "Second grouping"),
           keepdata = TRUE, keepmeta = TRUE)
}

check_pooled_gap <- function(x, family, details, test.overall) {
  forwarded <- captured <- component.labels <- NULL
  draw <- meta:::forest.meta
  testthat::local_mocked_bindings(
    forest.meta = function(...) {
      args <- list(...)
      # Hidden global summaries must not supply the component label width.
      args$text$common <- "C"
      args$text$random <- rep("R", length(args$x$method.random.ci))
      args$text$common.w <- rep("C", length(args$x$subgroup.levels))
      args$text$random.w <- rep("R", ncol(args$x$TE.random.w))
      if (family == "common") {
        component.labels <<- paste0(
          strrep("Long component common summary ", 4),
          args$x$subgroup.levels)
        args$text$common.w <- component.labels
      }
      else {
        component.labels <<- paste0(
          strrep("Long component random summary ", 4),
          colnames(args$x$TE.random.w))
        args$text$random.w <- component.labels
      }
      forwarded <<- args
      do.call(draw, args)
    },
    forest_meta_internal = function(...) {
      captured <<- list(...)
      invisible(NULL)
    }, .package = "meta")

  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off(), add = TRUE)
  forest(x, layout = "meta", common = TRUE, random = TRUE,
         prediction = FALSE, backtransf = FALSE, hetstat = FALSE,
         details = details, test.overall = test.overall,
         leftcols = c("studlab", "k", "effect", "ci"), rightcols = FALSE,
         calcwidth = list(common = TRUE, random = TRUE))

  expect_false(forwarded$overall)
  expect_true(all(is.na(forwarded$x$TE.common)))
  expect_true(all(is.na(forwarded$x$TE.random)))
  expect_true(all(forwarded$common.subgroup))
  expect_true(all(forwarded$random.subgroup))
  expect_identical(forwarded$details, details)
  expect_identical(forwarded$test$effect.subgroup, test.overall)
  expect_identical(captured$leftcols,
                   c("col.studlab", "col.k", "col.effect", "col.ci"))
  expect_false(captured$rsel)

  studlab <- captured$cols$col.studlab
  k <- captured$cols$col.k
  labels <- vapply(studlab$labels, function(g)
    if (is.character(g$label)) paste(g$label, collapse = "") else "", "")
  expected.count <- length(forwarded$x$subgroup.levels) *
    if (family == "common") 1L else length(forwarded$x$method.random.ci)
  component <- which(labels %in% component.labels & !is.na(studlab$rows))
  expect_length(component, expected.count)
  expect_setequal(labels[component], component.labels)

  k.left <- grid::convertWidth(sum(captured$x1[1:2]), "mm", valueOnly = TRUE)
  for (i in component) {
    # Numeric columns have compact label indices; match by rendered row.
    j <- match(studlab$rows[i], k$rows)
    expect_false(is.na(j), info = labels[i])
    expect_identical(trimws(as.character(k$labels[[j]]$label)), "8",
                     info = labels[i])
    label.width <- grid::convertWidth(grid::grobWidth(studlab$labels[[i]]),
                                      "mm", valueOnly = TRUE)
    expect_gte(k.left + 1e-7, label.width, label = labels[i])
  }
}

for (family in c("common", "random")) {
  for (multiple in c(FALSE, TRUE)) {
    for (details in c(FALSE, TRUE)) {
      for (test.overall in c(FALSE, TRUE)) {
        test_that(paste("component", family, "labels leave room for k with",
                        if (multiple) "multiple methods," else "one method,",
                        "details =", details, "and test.overall =", test.overall), {
          x <- pooled_gap_fixture(if (multiple) c("classic", "HK") else "classic")
          check_pooled_gap(x, family, details, test.overall)
        })
      }
    }
  }
}
