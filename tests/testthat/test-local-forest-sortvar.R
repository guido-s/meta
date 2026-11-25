test_that("sortvar data columns can share a name with a function", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- metacont(Ne, Me, Se, Nc, Mc, Sc, common = FALSE,
                data = Spooner2002, studlab = paste(author, year),
                subgroup = agegroup, print.subgroup.name = FALSE,
                label.e = "Nedocromil sodium", label.c = "Placebo")
  actual <- forest(x, sortvar = order)
  expected <- forest(x, sortvar = Spooner2002$order)
  expect_identical(actual$studlab, expected$studlab)
  expect_identical(actual$studlab,
                   x$studlab[order(x$subgroup, Spooner2002$order)])

  x <- metacont(Ne, Me, Se, Nc, Mc, Sc, data = Spooner2002,
                studlab = paste(author, year), subset = 1:10)
  actual <- forest(x, sortvar = order)
  expect_identical(actual$studlab,
                   x$studlab[order(Spooner2002$order[1:10])])
})

test_that("sortvar still resolves object fields, expressions, and caller vectors", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- metagen(c(0.3, 0.1, 0.2), rep(0.1, 3), studlab = LETTERS[1:3])
  expect_identical(forest(x)$studlab, x$studlab)
  expect_identical(forest(x, sortvar = TE)$studlab, x$studlab[order(x$TE)])
  expect_identical(forest(x, sortvar = -TE)$studlab, x$studlab[order(-x$TE)])
  ranks <- c(3, 2, 1)
  expect_identical(forest(x, sortvar = ranks)$studlab, rev(x$studlab))
  expect_error(forest(x, sortvar = 1:2), "have different length")
})

test_that("function-valued sortvar arguments give an informative error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (sel in list(TRUE, 1:10)) {
    x <- metacont(Ne, Me, Se, Nc, Mc, Sc, data = Spooner2002,
                  subset = sel)
    expect_error(forest(x, sortvar = sqrt),
                 "An object of class 'function' was supplied.",
                 fixed = TRUE)
  }
  dat <- Spooner2002
  dat$sqrt <- rev(seq_len(nrow(dat)))
  x <- metacont(Ne, Me, Se, Nc, Mc, Sc, data = dat,
                studlab = paste(author, year))
  expect_identical(forest(x, sortvar = sqrt)$studlab, rev(x$studlab))
})

test_that("sortvar rejects non-vector objects before length checks or subsetting", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  invalid <- list(data.frame(rank = 1:3), matrix(1:3), array(1:3, c(3, 1, 1)),
                  list(1, 2, 3), new.env())
  for (value in invalid) {
    x <- metagen(c(0.3, 0.1, 0.2), rep(0.1, 3))
    expect_error(forest(x, sortvar = value), "Argument 'sortvar' must be a vector")
    dat <- data.frame(TE = c(0.3, 0.1, 0.2), seTE = rep(0.1, 3))
    x <- metagen(TE, seTE, data = dat, subset = 1:2)
    x$data <- list(rank = value, .subset = c(TRUE, TRUE, FALSE))
    expect_error(forest(x, sortvar = rank), "Argument 'sortvar' must be a vector")
  }
})

test_that("sortvar accepts factors and dates", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- metagen(c(0.3, 0.1, 0.2), rep(0.1, 3), studlab = LETTERS[1:3])
  for (value in list(factor(c("c", "a", "b")),
                     as.Date("2026-01-01") + c(3, 1, 2)))
    expect_identical(forest(x, sortvar = value)$studlab, x$studlab[c(2, 3, 1)])
})
