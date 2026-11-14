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
