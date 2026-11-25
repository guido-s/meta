sortvar_fixture <- function(subset = TRUE) {
  dat <- data.frame(yi = c(.1, .4, .2, .7, .3, .9),
                     sei = c(.15, .12, .18, .14, .17, .13),
                     order = c(4, 1, 6, 2, 5, 3))
  metagen(yi, sei, data = dat, subset = subset,
           studlab = paste0("Study_", LETTERS[1:6]), method.tau = "DL")
}

sortvar_studies <- function(method, x, sortvar) {
  if (missing(sortvar))
    call <- substitute(FUN(x), list(FUN = as.name(method)))
  else
    call <- substitute(FUN(x, sortvar = KEY),
                       list(FUN = as.name(method), KEY = substitute(sortvar)))
  if (method == "print")
    x <- summary(x)
  output <- capture.output(result <- eval(call, list(x = x), parent.frame()))
  if (method == "print")
    return(sub(" .*", "", trimws(output[grepl("^Study_[A-F] +", output)])))
  labels <- result$studlab
  if (method == "estimates")
    labels <- labels[grepl("^Study_[A-F]$", labels)]
  sub("^.*(Study_[A-F]).*$", "\\1", labels)
}

for (method in c("metacum", "metainf", "estimates", "print")) {
  test_that(paste(method, "resolves sorting fields, expressions, and caller vectors"), {
    x <- sortvar_fixture()
    expect_identical(sortvar_studies(method, x), x$studlab)
    expect_identical(sortvar_studies(method, x, NULL), x$studlab)
    expect_identical(sortvar_studies(method, x, TE), x$studlab[order(x$TE)])
    expect_identical(sortvar_studies(method, x, -TE), x$studlab[order(-x$TE)])
    expect_identical(sortvar_studies(method, x, seTE), x$studlab[order(x$seTE)])
    expect_identical(sortvar_studies(method, x, order),
                     x$studlab[order(x$data$order)])
    ranks <- c(3, 6, 1, 5, 2, 4)
    expect_identical(sortvar_studies(method, x, ranks), x$studlab[order(ranks)])
    expect_identical(sortvar_studies(method, x, factor(ranks)),
                     x$studlab[order(ranks)])
    expect_identical(sortvar_studies(method, x, as.Date("2026-01-01") + ranks),
                     x$studlab[order(ranks)])
    expect_error(sortvar_studies(method, x, 1:2), "have different length")
  })

  test_that(paste(method, "sorts subsetted data and validates before subsetting"), {
    sel <- c(1, 3, 4, 6)
    x <- sortvar_fixture(sel)
    expect_identical(sortvar_studies(method, x, TE), x$studlab[order(x$TE)])
    expect_identical(sortvar_studies(method, x, order),
                     x$studlab[order(x$data$order[sel])])
    expect_identical(sortvar_studies(method, x, NULL), x$studlab)
    ranks <- c(3, 1, 4, 2)
    expect_identical(sortvar_studies(method, x, ranks), x$studlab[order(ranks)])
    for (value in list(matrix(6:1, ncol = 1), as.list(6:1))) {
      x$data$rank <- value
      expect_error(sortvar_studies(method, x, rank),
                   "Argument 'sortvar' must be a vector", fixed = TRUE)
    }
  })

  test_that(paste(method, "rejects unsupported sorting objects with informative errors"), {
    x <- sortvar_fixture()
    for (value in list(matrix(6:1, ncol = 1), array(6:1, c(6, 1, 1)),
                       as.list(6:1), data.frame(rank = 6:1), ~ TE, new.env()))
      expect_error(sortvar_studies(method, x, value),
                   "Argument 'sortvar' must be a vector", fixed = TRUE)
    expect_error(sortvar_studies(method, x, sqrt),
                 "An object of class 'function' was supplied.", fixed = TRUE)
  })
}
