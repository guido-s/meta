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
