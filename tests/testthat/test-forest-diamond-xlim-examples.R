test_that("forest clips pooled diamonds to narrow x-limits", {
  m <- metabin(1:3, 101:103, c(20, 10, 30), 101:103)
  xlim <- log(c(0.2, 5))
  assign(".forest_diamond_xs", list(), envir = .GlobalEnv)
  on.exit(rm(".forest_diamond_xs", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(grid::grid.polygon, quote({
      .forest_diamond_xs[[length(.forest_diamond_xs) + 1]] <<- as.numeric(x)
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::grid.polygon))),
          add = TRUE)
  #
  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = xlim))
  draw.ci.diamond(m$TE.common, m$lower.common, m$upper.common,
                  1, xlim[1], xlim[2], "gray", "black", 1)
  grid::popViewport()
  expect_silent(forest(m, xlim = c(0.2, 5)))
  dev.off()

  xs <- get(".forest_diamond_xs", envir = .GlobalEnv)
  expect_true(length(xs) > 0)
  expect_true(all(xs[[1]] >= xlim[1] & xs[[1]] <= xlim[2]))
})

test_that("forest accepts graphics devices", {
  data(Fleiss1993bin)
  m <- metabin(d.asp, n.asp, d.plac, n.plac, study,
               data = Fleiss1993bin, sm = "OR", method = "I")
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd), add = TRUE)
  on.exit(unlink(file.path(tempdir(), "Rplots.pdf")), add = TRUE)

  unlink("Rplots.pdf")
  expect_invisible(forest(m, device = pdf))
  expect_true(file.exists("Rplots.pdf"))

  unlink("Rplots.pdf")
  expect_invisible(forest(m, device = "pdf"))
  expect_true(file.exists("Rplots.pdf"))
})

test_that("BMJ layout centers combined treatment group headers", {
  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)

  fb <- forest(metabin(1:5, 101:105, 5:1, 101:105), layout = "BMJ")
  fc <- forest(metacont(1:5, rep(1, 5), 101:105,
                        5:1, rep(1, 5), 101:105),
               layout = "BMJ")
  fi <- forest(metainc(1:5, 101:105, 5:1, 101:105), layout = "BMJ")
  fin <- forest(metainc(1:5, 101:105, 5:1, 101:105,
                        n.e = 11:15, n.c = 21:25),
                layout = "BMJ")
  dev.off()

  expect_true(all(c("col.event.n.e", "col.event.n.c") %in% fb$leftcols))
  expect_true(all(c("col.mean.sd.n.e", "col.mean.sd.n.c") %in% fc$leftcols))
  expect_true(all(c("col.event.time.e", "col.event.time.c") %in% fi$leftcols))
  expect_true(all(c("col.event.time.n.e", "col.event.time.n.c") %in%
                    fin$leftcols))
})

test_that("forest accounts for attached treatment labels in column headings", {
  m <- metagen(1:5, 5:1, n.e = 1:5, n.c = 1:5, sm = "MD",
               label.e = "Experimental", label.c = "Control")

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  expect_silent(forest(m))
  expect_silent(
    forest(m, leftlabs = c(NA, NA, NA, "A\nTotal", NA))
  )
  dev.off()
})

test_that("forest accepts plot titles", {
  m <- metagen(1:5, 5:1, sm = "MD")

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  res <- forest(m, main = "First line\nSecond line", just.main = "left",
                xpos.main = 0, fs.main = 14, ff.main = "italic",
                col.main = "blue", gap.main = 2, lineheight.main = 1)

  expect_equal(res$main, "First line\nSecond line")
  expect_equal(res$just.main, "left")
  expect_equal(res$xpos.main, 0)
  expect_equal(res$fs.main, 14)
  expect_equal(res$ff.main, "italic")
  expect_equal(res$col.main, "blue")
  expect_equal(res$gap.main, 2)
  expect_equal(res$lineheight.main, 1)

  res.default <- forest(m, main = "Title", fontsize = 10)
  expect_equal(res.default$fs.main, 11)
  expect_error(forest(m, main = "Title", gap.main = 1.5),
               "Argument 'gap.main'")
  dev.off()
})

test_that("forest title does not shift header line", {
  m <- metagen(1:5, 5:1, sm = "MD")
  header.line.y <- character(0)
  assign(".forest_header_line_y", header.line.y, envir = .GlobalEnv)
  on.exit(rm(".forest_header_line_y", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(grid::grid.lines, quote({
      if (!missing(gp) && identical(gp$col, "hotpink"))
        .forest_header_line_y <<- c(.forest_header_line_y, as.character(y))
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::grid.lines))),
          add = TRUE)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  forest(m, header = TRUE, col.header.line = "hotpink")
  forest(m, main = "Title", header = TRUE, col.header.line = "hotpink")
  dev.off()

  header.line.y <- get(".forest_header_line_y", envir = .GlobalEnv)
  expect_true(length(header.line.y) > 0)
  expect_length(unique(header.line.y), 1)
})

test_that("forest accepts different colours for multiple overall diamonds", {
  m <- metagen(1:5, 5:1, sm = "MD")
  m$TE.common <- c(1, 2)
  m$lower.common <- c(0, 1)
  m$upper.common <- c(2, 3)
  m$statistic.common <- c(1, 2)
  m$pval.common <- c(0.1, 0.2)
  m$text.common <- c("A", "B")
  m$TE.random <- c(1.5, 2.5)
  m$lower.random <- c(0.5, 1.5)
  m$upper.random <- c(2.5, 3.5)
  m$statistic.random <- c(1.5, 2.5)
  m$pval.random <- c(0.15, 0.25)
  m$text.random <- c("C", "D")

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  expect_silent(
    forest(m, common = TRUE, random = TRUE,
           col.diamond.common = c("red", "blue"),
           col.diamond.random = c("green", "yellow"),
           col.diamond.lines.common = c("black", "gray"),
           col.diamond.lines.random = c("purple", "orange"))
  )
  expect_error(
    forest(m, common = TRUE, random = FALSE,
           col.diamond.common = c("red", "blue", "green")),
    "Length of argument 'col.diamond.common'"
  )
  dev.off()
})

test_that("forest accepts different colours for multiple prediction intervals", {
  m <- metagen(1:5, 5:1, method.predict = c("V", "HTS"))

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  expect_silent(
    forest(m, prediction = TRUE,
           col.predict = c("red", "blue"),
           col.predict.lines = c("green", "yellow"))
  )
  expect_error(
    forest(m, prediction = TRUE,
           col.predict = c("red", "blue", "green")),
    "Length of argument 'col.predict'"
  )
  dev.off()
})

test_that("forest accepts different colours for subgroup diamonds", {
  m <- metagen(1:5, 5:1, method.random.ci = c("classic", "HK"),
               subgroup = LETTERS[c(1, 1, 2, 2, 2)])

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  expect_silent(
    forest(m,
           col.diamond.common.subgroup = c("red", "blue"),
           col.diamond.random.subgroup = c("green", "yellow",
                                           "purple", "orange"),
           col.diamond.lines.common.subgroup = c("black", "gray"),
           col.diamond.lines.random.subgroup = c("blue", "green"))
  )
  expect_error(
    forest(m, col.diamond.random.subgroup = c("red", "blue", "green")),
    "Length of argument 'col.diamond.random.subgroup'"
  )
  dev.off()
})

test_that("subgroup diamond colours follow displayed subgroup order", {
  m <- metagen(1:5, 5:1, method.random.ci = c("classic", "HK"),
               subgroup = LETTERS[c(2, 2, 1, 1, 1)])

  expect_equal(
    setlength_subgroup_colors(
      1:4, 2, 4, sort(m$subgroup.levels), m$subgroup.levels,
      "number of random effects estimates", "col.diamond.random.subgroup"),
    c(3, 4, 1, 2)
  )
  expect_equal(
    setlength_subgroup_colors(
      1:4, 2, 4, m$subgroup.levels, m$subgroup.levels,
      "number of random effects estimates", "col.diamond.random.subgroup"),
    1:4
  )
})

test_that("forest accepts different colours for subgroup prediction intervals", {
  m <- metagen(1:5, 5:1, method.predict = c("V", "HTS"),
               subgroup = LETTERS[c(1, 1, 2, 2, 2)])

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  expect_silent(
    forest(m, prediction = TRUE, prediction.subgroup = TRUE,
           col.predict.subgroup = c("red", "blue"),
           col.predict.lines.subgroup = c("green", "yellow",
                                         "purple", "orange"))
  )
  expect_error(
    forest(m, prediction = TRUE, prediction.subgroup = TRUE,
           col.predict.subgroup = c("red", "blue", "green")),
    "Length of argument 'col.predict.subgroup'"
  )
  dev.off()
})

test_that("subgroup prediction interval colours follow displayed order", {
  m <- metagen(1:5, 5:1, method.predict = c("V", "HTS"),
               subgroup = LETTERS[c(2, 2, 1, 1, 1)])

  expect_equal(
    setlength_subgroup_colors(
      1:4, 2, 4, sort(m$subgroup.levels), m$subgroup.levels,
      "number of prediction intervals", "col.predict.subgroup"),
    c(3, 4, 1, 2)
  )
  expect_equal(
    setlength_subgroup_colors(
      1:4, 2, 4, m$subgroup.levels, m$subgroup.levels,
      "number of prediction intervals", "col.predict.subgroup"),
    1:4
  )
})

test_that("prediction.subgroup accepts method-specific subgroup matrices", {
  m <- metagen(1:5, 5:1, method.random.ci = c("classic", "HK"),
               subgroup = LETTERS[c(2, 2, 1, 1, 1)], prediction = TRUE,
               method.predict = c("V", "HK-PR", "S"))

  expect_equal(
    show_prediction_subgroup_results(
      m$df.predict.w > 1, length(m$subgroup.levels),
      m$lower.predict.w, m$upper.predict.w),
    c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE)
  )
  expect_equal(
    order_subgroup_results(
      show_prediction_subgroup_results(
        m$df.predict.w > 1, length(m$subgroup.levels),
        m$lower.predict.w, m$upper.predict.w),
      order(factor(m$subgroup.levels, levels = sort(m$subgroup.levels))),
      ncol(m$df.predict.w)),
    c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE)
  )

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  expect_silent(
    forest(m, col.diamond.random.subgroup = 1:4,
           prediction.subgroup = k.w > 2,
           col.predict = c("brown", "orange", "violet"))
  )
  expect_silent(
    forest(m, col.diamond.random.subgroup = 1:4,
           prediction.subgroup = df.predict.w > 1,
           col.predict = c("brown", "orange", "violet"))
  )
  dev.off()
})

test_that("random.subgroup accepts method-specific subgroup matrices", {
  m <- metagen(1:5, 5:1, method.random.ci = c("classic", "HK"),
               subgroup = LETTERS[c(2, 2, 1, 1, 1)], prediction = TRUE,
               method.predict = c("V", "HK-PR", "S"))

  expect_equal(
    show_subgroup_results_by_method(
      m$df.random.w > 1, length(m$subgroup.levels),
      m$lower.random.w, m$upper.random.w, "random.subgroup"),
    c(TRUE, FALSE, TRUE, TRUE)
  )
  expect_equal(
    order_subgroup_results(
      show_subgroup_results_by_method(
        m$df.random.w > 1, length(m$subgroup.levels),
        m$lower.random.w, m$upper.random.w, "random.subgroup"),
      order(factor(m$subgroup.levels, levels = sort(m$subgroup.levels))),
      ncol(m$df.random.w)),
    c(TRUE, TRUE, TRUE, FALSE)
  )

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  expect_silent(
    forest(m, random.subgroup = df.random.w > 1)
  )
  dev.off()
})
