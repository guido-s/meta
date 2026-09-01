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
