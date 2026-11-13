test_that("cidprop rounding lists preserve legacy settings and defaults", {
  m <- metagen(c(-0.4, -0.1, 0.2, 0.7), rep(0.2, 4), sm = "MD",
               method.tau = "DL", prediction = TRUE)
  x <- cidprop(m, cid = 0.25)
  before <- x
  expect_identical(capture.output(print(x)),
                   capture.output(print(x, digits = NULL)))
  expect_identical(
    capture.output(print(x, digits.cid = 3, digits.percent = 2)),
    capture.output(print(x, digits = list(cid = 3, percent = 2))))
  expect_identical(capture.output(print(x, digits.cid = 3)),
                   capture.output(print(x, digits = 3)))
  expect_identical(
    capture.output(print(x, digits.cid = 3, digits.percent = 2)),
    capture.output(print(x, digits = list(cid = 3), digits.percent = 2)))
  for (setting in c("cid", "percent")) {
    args <- setNames(list(2), paste0("digits.", setting))
    expect_error(do.call(print, c(list(x, digits = setNames(list(2), setting)),
                                  args)), "cannot both be specified")
  }
  expect_error(print(x, digits = 2, digits.cid = 2),
               "cannot both be specified")
  expect_error(print(x, digits = list(percent = -1)), "digits.percent")
  expect_error(print(x, digits = list(2)), "named list")
  expect_identical(x, before)
})

test_that("BLUP print lists preserve defaults, scalars and legacy arguments", {
  m <- metagen(c(-0.4, -0.1, 0.2, 0.7, -0.2, 0.3, 0.6, 0.9), rep(0.2, 8),
               sm = "RR", subgroup = rep(c("A", "B"), each = 4),
               method.tau = "DL", method.random.ci = "classic")
  x <- blup(m)
  before <- x
  expect_identical(capture.output(print(x)),
                   capture.output(print(x, print = NULL, digits = NULL)))
  for (backtransf in c(FALSE, TRUE)) {
    expect_identical(
      capture.output(print(x, backtransf = backtransf, se = TRUE,
                           digits = 3, digits.se = 2, digits.tau2 = 2,
                           digits.tau = 1, print.tau2 = TRUE, print.tau = TRUE)),
      capture.output(print(x, backtransf = backtransf, se = TRUE,
                           digits = list(digits = 3, se = 2, tau2 = 2, tau = 1),
                           print = list(tau2 = TRUE, tau = TRUE))))
  }
  expect_identical(capture.output(print(x, digits = 3)),
                   capture.output(print(x, digits = list(digits = 3))))
  expect_identical(
    capture.output(print(x, digits.se = 2, rowlab = paste0("Study ", 1:8))),
    capture.output(print(x, digits = list(se = 2),
                         rowlab = paste0("Study ", 1:8))))
  expect_identical(
    capture.output(print(x, digits.se = 2, digits.tau = 1,
                         print.tau = FALSE, print.tau2 = TRUE)),
    capture.output(print(x, digits = list(se = 2), digits.tau = 1,
                         print = list(tau = FALSE), print.tau2 = TRUE)))
  for (setting in c("se", "tau2", "tau")) {
    args <- setNames(list(2), paste0("digits.", setting))
    expect_error(do.call(print, c(list(x, digits = setNames(list(2), setting)),
                                  args)), "cannot both be specified")
  }
  for (setting in c("tau2", "tau")) {
    args <- setNames(list(TRUE), paste0("print.", setting))
    expect_error(do.call(print, c(list(x, print = setNames(list(TRUE), setting)),
                                  args)), "cannot both be specified")
  }
  expect_error(print(x, digits = list(se = -1)), "digits.se")
  expect_error(print(x, print = list(tau = NA)), "print.tau")
  expect_error(print(x, digits = list(2)), "named list")
  expect_identical(x, before)
})

test_that("estimates BLUP rounding defaults to object attributes", {
  m <- metagen(c(-0.4, -0.1, 0.2, 0.7), rep(0.2, 4), sm = "VE",
               method.tau = "DL", method.random.ci = "classic")
  for (backtransf in c(FALSE, TRUE)) {
    x <- estimates(blup(m), backtransf = backtransf, writexl = FALSE,
                   digits = 3, digits.se = 2, digits.tau2 = 1, digits.tau = 2)
    before <- x
    expect_identical(
      capture.output(print(x)),
      capture.output(print(x, digits = list(digits = 3, se = 2, tau2 = 1,
                                            tau = 2))))
    expect_identical(capture.output(print(x)),
                     capture.output(print(x, digits = NULL)))
    expect_identical(capture.output(print(x, digits = 4)),
                     capture.output(print(x, digits = list(digits = 4))))
    expect_identical(
      capture.output(print(x, digits.se = 3, rowlab = paste0("Study ", 1:4))),
      capture.output(print(x, digits = list(se = 3),
                           rowlab = paste0("Study ", 1:4))))
    changed <- x
    attr(changed, "digits") <- 4
    attr(changed, "digits.se") <- 3
    attr(changed, "digits.tau2") <- 2
    attr(changed, "digits.tau") <- 1
    expect_identical(
      capture.output(print(changed)),
      capture.output(print(x, digits = list(digits = 4, se = 3, tau2 = 2,
                                            tau = 1))))
    expect_identical(
      capture.output(print(changed)),
      capture.output(print(x, digits = 4, digits.se = 3, digits.tau2 = 2,
                           digits.tau = 1)))
    for (setting in c("se", "tau2", "tau")) {
      args <- setNames(list(2), paste0("digits.", setting))
      expect_error(do.call(print, c(list(x, digits = setNames(list(2), setting)),
                                    args)), "cannot both be specified")
    }
    expect_identical(x, before)
  }
  expect_error(print(x, digits = list(se = -1)), "digits.se")
  expect_error(print(x, digits = list(2)), "named list")
})

test_that("result print lists retain transformed and incomplete bound output", {
  for (sm in c("MD", "RR", "VE")) {
    m <- metagen(c(-0.4, -0.1, 0.2, 0.7), rep(0.2, 4), sm = sm,
                 method.tau = "DL", method.random.ci = "classic")
    x <- blup(m)
    x$lower[1] <- -Inf
    x$upper[1] <- Inf
    x$lower[2] <- NA_real_
    x$upper[3] <- NA_real_
    for (backtransf in c(FALSE, TRUE)) {
      expect_identical(
        capture.output(print(x, digits = 3, backtransf = backtransf)),
        capture.output(print(x, digits = list(digits = 3),
                             backtransf = backtransf)))
      y <- estimates(x, backtransf = backtransf, writexl = FALSE, digits = 3)
      expect_identical(capture.output(print(y)),
                       capture.output(print(y, digits = list(digits = 3))))
    }
    if (sm != "VE") {
      m$backtransf <- FALSE
      y <- cidprop(m, cid = if (sm == "RR") 0.8 else 0.25)
      expect_identical(capture.output(print(y, digits.cid = 3)),
                       capture.output(print(y, digits = list(cid = 3))))
    }
  }
})
