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
