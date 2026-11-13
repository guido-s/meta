test_that("long annotations reserve space without right-hand columns", {
  ma <- metagen(1:3, 3:1)
  plot_width <- function(...) {
    file <- tempfile(fileext = ".pdf")
    on.exit(unlink(file))
    forest(ma, rightcols = FALSE, filename = file, ...)$width
  }

  default.width <- plot_width()
  expect_gt(plot_width(xlab = strrep("Long label ", 12)), default.width)
  expect_gt(plot_width(label.right = strrep("Right label ", 12)),
            default.width)
})
