issue_forest_draw <- function(x, ...) {
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      observed <<- list(...)
      draw(...)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off())
  forest(x, ...)
  observed$widths.mm <- grid::convertWidth(observed$x1, "mm", valueOnly = TRUE)
  observed$label.widths.mm <- vapply(observed$cols$col.studlab$labels,
    function(g) grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
    numeric(1))
  observed
}

test_that("issue 108 defaults reserve real annotation width with custom columns", {
  x <- metaprop(c(13, 47, 15, 33, 35), c(129, 400, 103, 177, 38),
                method = "GLMM", common = FALSE, prediction = TRUE)
  for (layout in c("meta", "meta8", "RevMan5", "BMJ", "JAMA")) {
    for (totals in c(TRUE, FALSE)) {
      for (columns in list(NULL, "studlab", c("studlab", "n", "effect", "ci"),
                           c("n", "studlab", "effect", "ci"))) {
        args <- list(x = x, layout = layout, pooled.totals = totals,
                     print.I2.ci = TRUE, digits = 3)
        if (!is.null(columns)) {
          args$leftcols <- columns
          args$leftlabs <- ifelse(columns == "n", strrep("N", 50), columns)
          args$rightcols <- FALSE
        }
        drawn <- do.call(issue_forest_draw, args)
        labels <- lapply(drawn$cols$col.studlab$labels, `[[`, "label")
        het <- which(vapply(labels, function(z) is.call(z) &&
          is.character(z[[2]]) && startsWith(z[[2]], "Heterogeneity"), logical(1)))
        expect_length(het, 1)
        start <- 2 * match("col.studlab", drawn$leftcols) - 1
        available <- sum(drawn$widths.mm[seq.int(start, 2 * length(drawn$leftcols))])
        expect_gte(available + 1e-7, drawn$label.widths.mm[het],
                   label = paste(layout, totals, paste(columns, collapse = ",")))
      }
    }
  }
})

test_that("issue 109 native GLMM labels retain the selected distribution and full DF", {
  dat <- data.frame(event = c(10, 20, 30, 40, 50, 30, 40, 50, 60, 70),
                    total = 100, group = rep(c("A", "B"), each = 5))
  for (methods in list("HK", "classic", c("classic", "HK"), c("HK", "classic"))) {
    x <- metaprop(event, total, data = dat, sm = "PLOGIT", method = "GLMM",
                  method.random.ci = methods, common = FALSE,
                  subgroup = group, tau.common = TRUE)
    if ("HK" %in% methods) {
      i <- match("HK", methods)
      df <- if (length(methods) == 1) x$df.Q.b.random else x$df.Q.b.random[[i]]
      expect_equal(unname(df), c(1, 8))
      expect_equal(unname(x$pval.Q.b.random[i]),
                   unname(pf(x$Q.b.random[i], df[1], df[2], lower.tail = FALSE)))
    }
    for (layout in c("meta", "meta8", "RevMan5", "BMJ", "JAMA")) {
      drawn <- issue_forest_draw(
        x, layout = layout, label.test.subgroup.random = "ISSUE109",
        print.Q.subgroup = TRUE, digits.Q = 2, digits.pval.Q = 4,
        scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE)
      labels <- lapply(drawn$cols$col.studlab$labels, `[[`, "label")
      labels <- Filter(function(z) is.call(z) && identical(z[[2]], "ISSUE109"), labels)
      expect_length(labels, 1)
      label <- labels[[1]]
      hk <- methods[1] == "HK"
      if (layout %in% c("RevMan5", "BMJ")) {
        expected <- if (hk) "F" else if (layout == "BMJ") quote(chi^2) else quote("Chi"^2)
        expect_identical(label[[3]], expected)
        expect_match(label[[4]], if (hk) "df = 1, 8" else "df = 1$", fixed = hk)
        if (!hk)
          expect_false(grepl("df = 1,", label[[4]], fixed = TRUE))
      }
      else if (hk) {
        expect_identical(label[[3]][[2]], quote(F))
        expect_identical(label[[3]][[3]], "1, 8")
      }
      else {
        expect_identical(label[[3]], quote(chi[1L]^2))
      }
      if (hk) {
        text <- paste(deparse(label), collapse = " ")
        expect_match(text, "4.86", fixed = TRUE)
        expect_match(text, "0.0587", fixed = TRUE)
      }
    }
  }
})

test_that("custom column prefixes do not contribute to pooled label space", {
  x <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), n.e = c(10, 20, 30),
                method.tau = "DL", common = FALSE)
  for (layout in c("meta", "meta8", "RevMan5", "BMJ", "JAMA")) {
    for (totals in c(TRUE, FALSE)) {
      text <- strrep("Pooled label ", 15)
      drawn <- issue_forest_draw(
        x, layout = layout, leftcols = c("n.e", "studlab", "effect", "ci"),
        leftlabs = c(strrep("N", 50), "Study", "Effect", "CI"), rightcols = FALSE,
        pooled.totals = totals, text.random = text, calcwidth = list(random = TRUE),
        hetstat = FALSE, test.overall = FALSE, details = FALSE)
      labels <- lapply(drawn$cols$col.studlab$labels, `[[`, "label")
      pooled <- which(vapply(labels, identical, logical(1), text))
      expect_length(pooled, 1)
      # The effect column starts after the study column and its following gap.
      expect_gte(sum(drawn$widths.mm[3:4]) + 1e-7,
                 drawn$label.widths.mm[pooled], label = layout)
    }
  }
})

test_that("effect-test and heterogeneity options reserve only visible annotations", {
  x <- metagen(c(0.1, 0.3, 0.5, 0.2, 0.4, 0.6), rep(0.1, 6),
                subgroup = rep(c("A", "B"), each = 3), method.tau = "DL")
  tag <- strrep("EFFECT_TEST_", 20)
  for (layout in c("meta", "meta8", "RevMan5", "BMJ", "JAMA")) {
    for (enabled in c(TRUE, FALSE)) {
      drawn <- issue_forest_draw(
        x, layout = layout, leftcols = c("seTE", "studlab", "effect", "ci"),
        leftlabs = c(strrep("N", 50), "Study", "Effect", "CI"), rightcols = FALSE,
        test.overall = enabled, test.effect.subgroup = enabled,
        label.test.overall.common = tag, label.test.overall.random = tag,
        label.test.effect.subgroup.common = tag,
        label.test.effect.subgroup.random = tag, hetstat = enabled,
        print.I2.ci = TRUE, details = FALSE)
      labels <- lapply(drawn$cols$col.studlab$labels, `[[`, "label")
      selected <- which(vapply(labels, function(z) is.call(z) &&
        identical(z[[2]], tag), logical(1)))
      expect_length(selected, if (enabled) 6 else 0)
      if (enabled) {
        available <- sum(drawn$widths.mm[3:(2 * length(drawn$leftcols))])
        expect_gte(available + 1e-7, max(drawn$label.widths.mm[selected]), label = layout)
      }
    }
  }
})
