test_that("interval-slot identity does not affect effect sharing or printed columns", {
  x <- metabind_repeated_interval_fixture()
  for (model in c("random", "predict")) {
    rows <- x[[model]]$id == 2
    x[[model]]$interval.slot[rows] <- x[[model]]$interval.slot[rows] + 10L
    rows <- x$subgroup[[model]]$id == 2
    x$subgroup[[model]]$interval.slot[rows] <-
      x$subgroup[[model]]$interval.slot[rows] + 10L
  }
  prepared <- meta:::prepare_metabind_subgroups(x)
  expect_true(prepared$shared.overall)
  expect_true(prepared$same.results)
  expect_equal(nrow(prepared$subdat), sum(vapply(x$subgroup, nrow, integer(1))))
  expect_silent(output <- capture.output(print(x)))
  expect_false(any(grepl("interval.slot", output, fixed = TRUE)))
})

test_that("identical default labels retain repeated slots through preparation and mapping", {
  for (shared in c(FALSE, TRUE)) {
    x <- metabind_repeated_interval_fixture(labels = FALSE)
    if (!shared)
      x$samedata <- FALSE
    prepared <- meta:::prepare_metabind_subgroups(x)
    expect_identical(prepared$shared.overall, shared)
    dat <- prepared$dat
    for (type in c("Random effects model", "Prediction interval"))
      expect_equal(sum(dat$type == type), if (shared) 2L else 4L)
    names(dat) <- sub("\\.pooled$", "", names(dat))
    skeleton <- metagen(rep(0, nrow(prepared$subdat)), rep(1, nrow(prepared$subdat)),
                         subgroup = as.character(prepared$subdat$id),
                         sm = x$overall$sm[1], common = FALSE, random = FALSE,
                         method.tau = "DL", method.tau.ci = "", warn = FALSE)
    mapped <- meta:::map_metabind_subgroups(
      skeleton, x, dat, dat, shared.overall = shared,
      common = TRUE, random = TRUE, prediction = TRUE, overall = TRUE,
      shared.heterogeneity = prepared$shared.heterogeneity)
    expect_length(mapped$ma$method.random.ci, 2L)
    expect_length(mapped$ma$method.predict, 2L)
    expect_equal(ncol(mapped$ma$TE.random.w), 2L)
    expect_equal(ncol(mapped$ma$lower.predict.w), 2L)
  }
})

test_that("different occurrence counts do not share an incomplete result set", {
  x <- metabind_repeated_interval_fixture(labels = FALSE)
  for (model in c("random", "predict")) {
    x[[model]] <- subset(x[[model]], !(id == 2 & interval.slot == 2))
    x$subgroup[[model]] <- subset(x$subgroup[[model]], !(id == 2 & interval.slot == 2))
  }
  prepared <- meta:::prepare_metabind_subgroups(x)
  expect_false(prepared$shared.overall)
  expect_false(prepared$same.results)
  expect_equal(sum(prepared$dat$type == "Random effects model"), 3L)
  expect_equal(sum(prepared$dat$type == "Prediction interval"), 3L)
})
