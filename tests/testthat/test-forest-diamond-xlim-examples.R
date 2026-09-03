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

test_that("settings.meta rejects NULL for non-null numeric settings", {
  on.exit(settings.meta(reset = TRUE), add = TRUE)

  expect_error(settings.meta(fontsize = NULL),
               "Argument 'fontsize' must not be NULL.", fixed = TRUE)
  expect_silent(settings.meta(fs.hetstat = NULL))
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

test_that("forest applies colours to individual study-row text", {
  m <- metagen(1:3, 1:3, sm = "MD")
  m$exclude <- c(FALSE, TRUE, FALSE)
  old.settings <- settings.meta(quietly = TRUE)
  on.exit(settings.meta(old.settings), add = TRUE)
  assign(".forest_study_text_colours", character(), envir = .GlobalEnv)
  on.exit(rm(".forest_study_text_colours", envir = .GlobalEnv), add = TRUE)
  suppressMessages(capture.output(
    trace(grid::grid.draw, quote({
      if (inherits(x, "text") && !is.null(x$gp$col))
        .forest_study_text_colours <<-
          c(.forest_study_text_colours, as.character(x$gp$col))
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::grid.draw))),
          add = TRUE)

  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f), add = TRUE)

  forest(m, filename = f, col.study.text = c("red", "blue", "green"))
  settings.meta(col.study.text = c("red", "blue", "green"), quietly = TRUE)
  forest(m, filename = f)

  colours <- get(".forest_study_text_colours", envir = .GlobalEnv)
  expect_true(all(c("red", "blue", "green", "lightgrey") %in% colours))
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

test_that("large file-output titles reserve enough rows", {
  m <- metagen(1:5, 5:1, sm = "MD")
  file1 <- tempfile(fileext = ".pdf")
  file2 <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(file1, file2)), add = TRUE)

  res1 <- forest(m, main = "MY\nTITLE", filename = file1)
  res2 <- forest(m, main = "MY\nTITLE", fs.main = 30,
                 header.line = "both", filename = file2)

  expect_true(file.exists(file2))
  expect_gt(res2$height, res1$height)
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

test_that("forest prints non-study labels on the left", {
  assign(".forest_add_text", list(), envir = .GlobalEnv)
  on.exit(rm(".forest_add_text", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(meta:::add.text, quote({
      .forest_add_text[[length(.forest_add_text) + 1]] <<-
        vapply(x$labels,
               function(z) paste(as.character(z$label), collapse = ""),
               character(1))
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(meta:::add.text))),
          add = TRUE)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  res <- forest(metagen(1:3, 1:3),
                leftcols = c("TE", "seTE"),
                rightcols = c("studlab", "effect", "ci"),
                details = TRUE)
  res2 <- forest(metagen(1:3, 1:3), leftcols = c("TE", "seTE"),
                 leftlabs = c("a", "b"), hetstat = FALSE)
  dev.off()

  drawn <- get(".forest_add_text", envir = .GlobalEnv)
  labels <- unlist(drawn)
  right.studlab <- drawn[[which(vapply(drawn, function(x) "Study" %in% x,
                                      logical(1)))]]

  expect_equal(res$leftcols, c("col.TE", "col.seTE"))
  expect_equal(res$rightcols, c("col.studlab", "col.effect", "col.ci"))
  expect_true(any(grepl("Heterogeneity", labels)))
  expect_true(any(grepl("Details of meta-analysis methods", labels)))
  expect_false(any(grepl("Heterogeneity", right.studlab)))
  expect_false(any(grepl("Details of meta-analysis methods", right.studlab)))

  expect_equal(res2$leftcols, c("col.TE", "col.seTE"))
  expect_equal(as.character(res2$colgap.studlab), "0mm")
})

test_that("forest suppresses non-study labels if no left columns are shown", {
  assign(".forest_add_text", list(), envir = .GlobalEnv)
  on.exit(rm(".forest_add_text", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(meta:::add.text, quote({
      .forest_add_text[[length(.forest_add_text) + 1]] <<-
        vapply(x$labels,
               function(z) paste(as.character(z$label), collapse = ""),
               character(1))
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(meta:::add.text))),
          add = TRUE)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  res <- forest(metagen(1:3, 1:3), leftcols = FALSE,
                rightcols = c("studlab", "effect", "ci"), details = TRUE)
  dev.off()

  labels <- unlist(get(".forest_add_text", envir = .GlobalEnv))

  expect_equal(res$leftcols, character(0))
  expect_equal(res$rightcols, c("col.studlab", "col.effect", "col.ci"))
  expect_false(any(grepl("Heterogeneity", labels)))
  expect_false(any(grepl("Details of meta-analysis methods", labels)))
})

test_that("calcwidth.pooled reserves extra spacing for pooled labels", {
  assign(".forest_widths", list(), envir = .GlobalEnv)
  on.exit(rm(".forest_widths", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(grid::grid.layout, quote({
      if (!missing(widths))
        .forest_widths[[length(.forest_widths) + 1]] <<- as.character(widths)
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::grid.layout))),
          add = TRUE)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  forest(metagen(1:3, 1:3, random = FALSE),
         leftcols = c("TE", "seTE"), leftlabs = c("b", "c"),
         hetstat = FALSE, calcwidth.pooled = FALSE)
  forest(metagen(1:3, 1:3, random = FALSE),
         leftcols = c("TE", "seTE"), leftlabs = c("b", "c"),
         hetstat = FALSE, calcwidth.pooled = TRUE)
  forest(metacont(101:103, 1:3, rep(1, 3),
                  201:203, 6:4, rep(1, 3)),
         leftcols = c("studlab", "n.e", "n.c"), calcwidth.pooled = FALSE)
  dat <- data.frame(TE = 1:3, seTE = 1:3, grp = c("a", "b", "c"))
  ma <- metagen(TE, seTE, data = dat, random = FALSE)
  forest(ma, leftcols = c("studlab", "grp"), leftlabs = c("", "Group"),
         hetstat = FALSE, calcwidth.pooled = TRUE)
  forest(ma, leftcols = c("studlab", "grp"), leftlabs = c("", "Group"),
         data.pooled = data.frame(grp = "Common effect model"),
         hetstat = FALSE, calcwidth.pooled = TRUE)
  dev.off()

  widths <- get(".forest_widths", envir = .GlobalEnv)
  expect_length(widths, 5)
  expect_false(any(grepl("sum(0mm, max(0mm", widths[[1]], fixed = TRUE)))
  expect_true(any(grepl("sum(0mm, max(0mm", widths[[2]], fixed = TRUE)))
  expect_equal(lengths(gregexpr("grobwidth", widths[[3]][1], fixed = TRUE)),
               4)
  expect_equal(lengths(gregexpr("-1*max", widths[[4]][2], fixed = TRUE)),
               2)
  expect_equal(lengths(gregexpr("-1*max", widths[[5]][2], fixed = TRUE)),
               1)

  ma <- metagen(1:5, 1:5)
  file1 <- tempfile(fileext = ".pdf")
  file2 <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(file1, file2)), add = TRUE)
  args <- list(header = "bo", details = TRUE, main = "MY TITLE",
               xlab = "djdjdjdjjd", label.left = "djdjdj",
               label.right = "jdjdjdh", layout = "R",
               overall.hetstat = FALSE)
  res1 <- do.call(forest, c(list(x = ma, filename = file1,
                                 calcwidth.pooled = FALSE), args))
  res2 <- do.call(forest, c(list(x = ma, filename = file2,
                                  calcwidth.pooled = TRUE), args))
  expect_gt(res2$width, res1$width)

  data(Olkin1995)
  meta1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
                   data = Olkin1995, subset = c(41, 47, 51, 59),
                   sm = "RR", method = "I",
                   studlab = paste(author, year))
  file3 <- tempfile(fileext = ".pdf")
  file4 <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(file3, file4)), add = TRUE)
  args <- list(x = meta1,
               leftcols = c("studlab", "n.e", "event.e", "n.c", "event.c"))
  res3 <- do.call(forest, c(list(filename = file3,
                                 calcwidth.pooled = FALSE), args))
  res4 <- do.call(forest, c(list(filename = file4,
                                 calcwidth.pooled = TRUE), args))
  expect_gt(res4$width, res3$width)

  meta2 <- update(meta1,
                  subgroup = ifelse(year < 1987,
                                    "Before 1987", "1987 and later"),
                  print.subgroup.name = FALSE)
  assign(".forest_widths", list(), envir = .GlobalEnv)
  file5 <- tempfile(fileext = ".pdf")
  on.exit(unlink(file5), add = TRUE)
  forest(meta2, filename = file5, colgap.studlab = "0mm",
         calcwidth.tests = TRUE)

  widths <- get(".forest_widths", envir = .GlobalEnv)
  expect_gt(length(widths), 0)
  expect_true(grepl("sum(0mm", widths[[1]][2], fixed = TRUE))
  expect_equal(widths[[1]][10], "2mm")
})

test_that("forest accepts multi-line labels", {
  assign(".forest_text", character(0), envir = .GlobalEnv)
  on.exit(rm(".forest_text", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(grid::grid.text, quote({
      .forest_text <<- c(.forest_text, as.character(label))
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::grid.text))),
          add = TRUE)
  suppressMessages(capture.output(
    trace(grid::grid.draw, quote({
      if (inherits(x, "text"))
        .forest_text <<- c(.forest_text, as.character(x$label))
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::grid.draw))),
          add = TRUE)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  forest(metagen(1:3, 1:3),
         smlab = "Summary 1\nSummary 2\nSummary 3",
         xlab = "X 1\nX 2\nX 3",
         label.left = "Left 1\nLeft 2\nLeft 3",
         label.right = "Right 1\nRight 2\nRight 3",
         bottom.lr = TRUE)
  forest(metagen(1:3, 1:3),
         label.left = "Top left 1\nTop left 2\nTop left 3",
         label.right = "Top right 1\nTop right 2\nTop right 3")
  dev.off()

  text <- get(".forest_text", envir = .GlobalEnv)
  expect_true(all(c("Summary 1", "Summary 2", "Summary 3",
                    "X 1", "X 2", "X 3",
                    "Left 1", "Left 2", "Left 3",
                    "Right 1", "Right 2", "Right 3",
                    "Top left 1", "Top left 2", "Top left 3",
                    "Top right 1", "Top right 2", "Top right 3") %in% text))
})

test_that("multi-line top labels reserve header rows", {
  assign(".forest_rows", data.frame(label = character(), row = integer()),
         envir = .GlobalEnv)
  assign(".forest_row", NA_integer_, envir = .GlobalEnv)
  on.exit(rm(".forest_rows", ".forest_row", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(grid::pushViewport, quote({
      vp <- list(...)[[1]]
      if (inherits(vp, "viewport") && !is.null(vp$layout.pos.row))
        .forest_row <<- vp$layout.pos.row[1]
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::pushViewport))),
          add = TRUE)
  suppressMessages(capture.output(
    trace(grid::grid.draw, quote({
      if (inherits(x, "text"))
        .forest_rows <<-
          rbind(.forest_rows,
                data.frame(label = as.character(x$label), row = .forest_row))
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::grid.draw))),
          add = TRUE)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  forest(metagen(1:3, 1:3),
         smlab = "p\nq\nr\ns\nt",
         main = "MY TITLE !!!\n:-)")
  dev.off()

  rows <- get(".forest_rows", envir = .GlobalEnv)
  expect_equal(rows$row[match(letters[16:20], rows$label)], 4:8)
  expect_gt(min(rows$row[rows$label %in% c("1.0000", "2.0000", "3.0000")]),
            8)

  assign(".forest_rows", data.frame(label = character(), row = integer()),
         envir = .GlobalEnv)
  assign(".forest_row", NA_integer_, envir = .GlobalEnv)

  pdf(tempfile())
  forest(metagen(1:3, 1:3),
         label.left = "M\nN\nO", label.right = "P\nQ\nR",
         bottom.lr = FALSE, header.line = "both")
  dev.off()

  rows <- get(".forest_rows", envir = .GlobalEnv)
  expect_equal(rows$row[match(LETTERS[13:15], rows$label)], 1:3)
  expect_equal(rows$row[match(LETTERS[16:18], rows$label)], 1:3)
  expect_gt(min(rows$row[rows$label %in% c("1.0000", "2.0000", "3.0000")]),
            3)
})

test_that("bottom annotations use structured rows", {
  assign(".axis_rows", list(), envir = .GlobalEnv)
  assign(".text_rows", data.frame(label = character(), row = integer()),
         envir = .GlobalEnv)
  on.exit(rm(".axis_rows", ".text_rows", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(meta:::draw.axis, quote({
      .axis_rows[[length(.axis_rows) + 1]] <<-
        c(axis = axis.row, labels = axis.label.row)
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(meta:::draw.axis))),
          add = TRUE)
  suppressMessages(capture.output(
    trace(meta:::add.text, quote({
      labs <- vapply(x$labels,
                     function(z) paste(as.character(z$label), collapse = ""),
                     character(1))
      keep <- grepl("Heterogeneity|Test for overall|Du$|^a$", labs)
      if (any(keep))
        .text_rows <<-
          rbind(.text_rows,
                data.frame(label = labs[keep], row = x$rows[keep]))
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(meta:::add.text))),
          add = TRUE)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  forest(metagen(1:5, 1:5), xlab = "a", label.left = "Du",
         test.overall = TRUE)
  dev.off()

  axis.rows <- get(".axis_rows", envir = .GlobalEnv)[[1]]
  text.rows <- get(".text_rows", envir = .GlobalEnv)
  het.row <- text.rows$row[grep("Heterogeneity", text.rows$label)[1]]
  test.rows <- text.rows$row[grep("Test for overall", text.rows$label)]
  label.row <- text.rows$row[text.rows$label == "Du"][1]
  xlab.row <- text.rows$row[text.rows$label == "a"][1]

  expect_equal(unname(axis.rows["labels"]), het.row)
  expect_equal(test.rows, het.row + 1:2)
  expect_equal(label.row, unname(axis.rows["labels"]) + 1)
  expect_equal(xlab.row, label.row + 1)
})

test_that("heterogeneity statistics and axis use smaller default font size", {
  assign(".hetstat_fs", numeric(), envir = .GlobalEnv)
  assign(".test_fs", numeric(), envir = .GlobalEnv)
  assign(".addline_fs", numeric(), envir = .GlobalEnv)
  assign(".details_fs", numeric(), envir = .GlobalEnv)
  assign(".axis_fs", numeric(), envir = .GlobalEnv)
  on.exit(rm(".hetstat_fs", ".test_fs", ".addline_fs", ".details_fs",
             ".axis_fs", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(meta:::tg, quote({
      lab <- paste(as.character(x), collapse = "")
      if (grepl("Heterogeneity", lab))
        .hetstat_fs <<- c(.hetstat_fs, fs)
      if (grepl("Test for overall", lab))
        .test_fs <<- c(.test_fs, fs)
      if (lab == "ADD")
        .addline_fs <<- c(.addline_fs, fs)
      if (grepl("^- ", lab))
        .details_fs <<- c(.details_fs, fs)
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(meta:::tg))), add = TRUE)
  suppressMessages(capture.output(
    trace(grid::grid.text, quote({
      if (!missing(label) && is.numeric(label) && length(label) > 1 &&
          !missing(gp))
        .axis_fs <<- c(.axis_fs, gp$fontsize)
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::grid.text))),
          add = TRUE)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  forest(metagen(1:5, 1:5), fontsize = 10, test.overall = TRUE,
         text.addline1 = "ADD", details = TRUE)
  dev.off()

  expect_true(any(get(".hetstat_fs", envir = .GlobalEnv) == 9))
  expect_true(any(get(".test_fs", envir = .GlobalEnv) == 9))
  expect_true(any(get(".addline_fs", envir = .GlobalEnv) == 9))
  expect_true(any(get(".details_fs", envir = .GlobalEnv) == 9))
  expect_equal(get(".axis_fs", envir = .GlobalEnv)[[1]], 9)
})

test_that("forest draws short explicit x-axis ticks", {
  assign(".tick_lengths", numeric(), envir = .GlobalEnv)
  on.exit(rm(".tick_lengths", envir = .GlobalEnv), add = TRUE)

  suppressMessages(capture.output(
    trace(grid::grid.segments, quote({
      if (!missing(y0) && !missing(y1))
        .tick_lengths <<-
          c(.tick_lengths,
            abs(as.numeric(grid::convertY(y1 - y0, "lines",
                                          valueOnly = TRUE))))
    }), print = FALSE)
  ))
  on.exit(suppressMessages(capture.output(untrace(grid::grid.segments))),
          add = TRUE)

  pdf(tempfile())
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  forest(metagen(1:5, 1:5), xlab = "a", label.left = "Du")
  dev.off()

  expect_true(length(get(".tick_lengths", envir = .GlobalEnv)) > 0)
  expect_true(all(get(".tick_lengths", envir = .GlobalEnv) == 0.25))
})

test_that("file output handles forest layout height", {
  save_forest <- function(x, ...) {
    file <- tempfile(fileext = ".pdf")
    on.exit(unlink(file), add = TRUE)
    forest(x, ..., filename = file)
  }

  res1 <- save_forest(metagen(1:3, 1:3))
  res2 <- save_forest(metagen(1:3, 1:3),
                      xlab = "a\nb\nc", label.left = "A\nB\nC")
  res3 <- save_forest(metagen(1:3, 1:3), xlab = "a\nb\nc")
  ma <- metagen(1:5, 1:5)
  res4 <- save_forest(ma, details = TRUE)
  res5 <- save_forest(ma, details = TRUE, xlab = "A")
  res6 <- save_forest(ma, detail = TRUE, main = "ABC", fs.main = 30)
  res7 <- save_forest(ma, main = "MY\nTITLE", fs.main = 26,
                      detail = TRUE, hetstat = FALSE)
  res10 <- save_forest(ma, main = "MY\nTITLE", fs.main = 26,
                       detail = TRUE)
  res11 <- save_forest(ma, main = "TITLE")
  res12 <- save_forest(ma, main = "TITLE", fs.axis = 20)
  res13 <- save_forest(ma, main = "TITLE", xlab = "A", fs.axis = 20)
  res14 <- save_forest(ma, main = "TITLE", xlab = "A", fs.axis = 20,
                       fs.xlab = 45)

  expect_gt(res2$height, res1$height)
  expect_gte(res3$height, res1$height)
  expect_equal(res5$height, res4$height, tolerance = 0.001)
  expect_gt(res6$height, res4$height)
  expect_lt(res6$height - res4$height, 1)
  expect_gt(res10$height, res7$height)
  expect_gt(res11$height, save_forest(ma)$height)
  expect_gt(res12$height, res11$height)
  expect_gt(res13$height, save_forest(ma, main = "TITLE", xlab = "A")$height)
  expect_gt(res14$height, res13$height)
})

test_that("forest_dims matches file output dimensions", {
  ma <- metagen(1:5, 1:5)
  cases <- list(
    list(),
    list(xlab = "a\nb\nc"),
    list(details = TRUE),
    list(details = TRUE, xlab = "A"),
    list(main = "ABC", fs.main = 30)
  )

  for (args in cases) {
    file <- tempfile(fileext = ".pdf")
    on.exit(unlink(file), add = TRUE)
    dims <- do.call(forest_dims, c(list(x = ma), args))
    res <- do.call(forest, c(list(x = ma, filename = file), args))

    expect_equal(dims$width, res$width, tolerance = 0.001)
    expect_equal(dims$height, res$height, tolerance = 0.001)
  }
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
