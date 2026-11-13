interval_slot_fixture <- function(labels = FALSE, shared = TRUE) {
  ma <- metagen(c(-0.4, 0.1, 0.7, 0.3, 0.9, 1.2, 0.2, 0.6),
                 rep(c(0.15, 0.2), 4), sm = "MD",
                 subgroup = rep(c("A", "B"), each = 4),
                 method.tau = "DL", method.random.ci = c("HK", "HK"),
                 adhoc.hakn.ci = "se", prediction = TRUE,
                 method.predict = c("V", "V"))
  if (labels) {
    ma$text.random <- c("First CI", "Second CI")
    ma$text.predict <- c("First PI", "Second PI")
  }
  x <- metabind(ma, ma, keepmeta = TRUE)
  if (!shared)
    x$data <- list(x$data, x$data)
  x
}

interval_slot_draw <- function(x, ...) {
  observed <- NULL
  draw <- meta:::forest.meta
  testthat::local_mocked_bindings(forest.meta = function(...) {
    observed <<- list(...)
    do.call(draw, observed)
  }, .package = "meta")
  grDevices::pdf(NULL, width = 24, height = 18)
  on.exit(grDevices::dev.off())
  result <- forest(x, layout = "meta", ...)
  list(args = observed, result = result)
}

test_that("real rendering retains repeated default and custom-labelled intervals", {
  for (labels in c(FALSE, TRUE)) for (shared in c(FALSE, TRUE)) {
    x <- interval_slot_fixture(labels, shared)
    # Same method descriptors must not cause the second payload to overwrite one.
    for (family in c("random", "predict")) {
      second <- x[[family]]$interval.slot == 2L
      x[[family]]$lower[second] <- x[[family]]$lower[second] - 0.25
      x[[family]]$upper[second] <- x[[family]]$upper[second] + 0.25
    }
    before <- x
    expect_silent(output <- capture.output(print(x, overall = TRUE)))
    expect_false(any(grepl("interval.slot|interval.occurrence", output)))
    expect_silent(drawn <- interval_slot_draw(x, backtransf = FALSE))
    native <- drawn$args$x
    for (family in c("random", "predict")) {
      expect_length(native[[paste0("method.", if (family == "random")
        "random.ci" else "predict")]], 2L)
      for (id in if (shared) 1L else 1:2) {
        rows <- x[[family]]$id == id
        for (field in c("lower", "upper")) {
          actual <- if (shared) native[[paste0(field, ".", family)]] else
            unname(native[[paste0(field, ".", family, ".w")]][id, ])
          expect_equal(actual, x[[family]][[field]][rows])
        }
      }
    }
    expect_identical(x, before)
  }
})

test_that("slot joins reject ambiguous parents and survive independent row reversals", {
  x <- interval_slot_fixture(TRUE)
  for (family in c("random", "predict")) {
    bad <- x
    bad[[family]] <- rbind(bad[[family]], bad[[family]][1, ])
    expect_error(as.data.frame(bad), "Interval slots must uniquely identify")
    x[[family]] <- x[[family]][rev(seq_len(nrow(x[[family]]))), ]
    x$subgroup[[family]] <- x$subgroup[[family]][
      rev(seq_len(nrow(x$subgroup[[family]]))), ]
  }
  dat <- as.data.frame(x, drop_all_NA = FALSE)
  expect_equal(nrow(dat), sum(vapply(x$subgroup, nrow, 0L)))
  for (family in c("random", "predict")) {
    rows <- dat$type == if (family == "random")
      "Random effects model" else "Prediction interval"
    parent <- x[[family]]
    index <- match(paste(dat$id[rows], dat$interval.slot[rows]),
                   paste(parent$id, parent$interval.slot))
    expect_equal(dat$lower[rows], parent$lower[index])
    expect_equal(dat$text[rows], parent$text[index])
  }
  expect_true(meta:::prepare_metabind_subgroups(x)$shared.overall)
  expect_silent(drawn <- interval_slot_draw(x, backtransf = FALSE))
  expect_equal(drawn$args$x$lower.random, x$random$lower[1:2])
})

test_that("method ordering and slot renumbering do not change sharing", {
  x <- interval_slot_fixture(TRUE)
  for (family in c("random", "predict")) {
    rows <- x[[family]]$id == 2L
    x[[family]]$interval.slot[rows] <- rev(x[[family]]$interval.slot[rows]) + 10L
    rows <- x$subgroup[[family]]$id == 2L
    x$subgroup[[family]]$interval.slot[rows] <-
      13L - x$subgroup[[family]]$interval.slot[rows]
  }
  expect_true(meta:::prepare_metabind_subgroups(x)$shared.overall)
  for (family in c("random", "predict")) {
    selected <- x
    selected$common <- NULL
    other <- if (family == "random") "predict" else "random"
    selected[[other]] <- NULL
    for (part in c("pooled", "subgroup")) {
      dat <- if (part == "pooled") selected[[family]] else selected$subgroup[[family]]
      dat <- dat[dat$interval.slot %in% c(2L, 11L), ]
      if (part == "pooled") selected[[family]] <- dat else
        selected$subgroup[[family]] <- dat
    }
    expect_true(meta:::prepare_metabind_subgroups(selected)$shared.overall)
    expect_silent(drawn <- interval_slot_draw(selected, backtransf = FALSE))
    expect_length(drawn$args$x[[paste0("lower.", family)]], 1L)
    expect_silent(capture.output(print(selected)))
  }
})

test_that("vector DF and variance components survive print and real forest", {
  x <- interval_slot_fixture(shared = FALSE)
  fields <- c("tau2", "lower.tau2", "upper.tau2", "se.tau2",
              "tau", "lower.tau", "upper.tau")
  for (field in fields)
    x$overall[[field]] <- list(c(between = 0.1, within = 0.2),
                               c(between = 0.3, within = 0.4))
  x$overall$three.level <- TRUE
  x$subgroup$random$df.Q.b <- rep(list(c(1, 7.5), c(1, 9.25)), 4)
  x$subgroup$random$test.Q.b <- "F"
  before <- x
  raw <- as.data.frame(x, drop_all_NA = FALSE)
  for (field in fields)
    expect_identical(raw[[field]], x$overall[[field]][raw$id])
  expect_silent(output <- capture.output(print(
    x, overall = TRUE, print = list(tau2 = TRUE, tau = TRUE))))
  expect_true(any(grepl("1, 7.50", output, fixed = TRUE)))
  ungrouped <- x
  ungrouped$subgroup <- NULL
  expect_silent(variance.output <- capture.output(print(
    ungrouped, print = list(tau2 = TRUE, tau = TRUE),
    digits = list(tau2 = 4, tau = 4))))
  expect_true(any(grepl("0.1000, 0.2000", variance.output, fixed = TRUE)))
  expect_silent(drawn <- interval_slot_draw(
    x, print = list(tau2 = TRUE, tau = TRUE, tau2.ci = TRUE, tau.ci = TRUE)))
  expect_equal(unname(drawn$args$x$tau2.w), c(0.3, 0.7))
  expect_equal(unname(drawn$args$x$tau.w), sqrt(c(0.05, 0.25)))
  expect_true(all(is.na(drawn$args$x$lower.tau2.w)))
  expect_true(any(drawn$result$data$df.Q.b == "1, 7.50"))
  expect_identical(x, before)
})

test_that("occurrence ranks handle missing IDs and preserve within-family multiplicity", {
  dat <- data.frame(id = c(NA_integer_, NA_integer_, 2L, 2L),
                    interval.slot = c(2L, 1L, 12L, 11L))
  descriptors <- data.frame(method = rep("HK", 4))
  expect_identical(meta:::metabind_interval_occurrence(dat, descriptors),
                   c(2L, 1L, 2L, 1L))
})

test_that("identical repeated payloads render twice without exposing slot metadata", {
  x <- interval_slot_fixture()
  expect_silent(drawn <- interval_slot_draw(x))
  expect_length(drawn$args$x$lower.random, 2L)
  expect_length(drawn$args$x$lower.predict, 2L)
  expect_equal(drawn$args$x$lower.random, x$random$lower[1:2])
  expect_equal(drawn$args$x$lower.predict, x$predict$lower[1:2])
  expect_false(any(grepl("interval.slot|interval.occurrence",
                         c(drawn$result$leftcols, drawn$result$rightcols,
                           drawn$args$x$.text.details.methods))))
})

test_that("different source method orders map by descriptors rather than slot numbers", {
  ma <- interval_slot_fixture()$meta[[1]]
  first <- update(ma, method.random.ci = c("classic", "HK"),
                   method.predict = c("V", "HK"))
  second <- update(ma, method.random.ci = c("HK", "classic"),
                    method.predict = c("HK", "V"))
  for (shared in c(FALSE, TRUE)) {
    x <- metabind(first, second)
    if (!shared)
      x$data <- list(x$data, x$data)
    for (family in c("random", "predict")) {
      rows <- x[[family]]$id == 2L
      x[[family]]$interval.slot[rows] <- x[[family]]$interval.slot[rows] + 100L
      rows <- x$subgroup[[family]]$id == 2L
      x$subgroup[[family]]$interval.slot[rows] <-
        x$subgroup[[family]]$interval.slot[rows] + 100L
    }
    expect_identical(meta:::prepare_metabind_subgroups(x)$shared.overall, shared)
    expect_silent(drawn <- interval_slot_draw(x, backtransf = FALSE))
    for (family in c("random", "predict")) {
      method <- drawn$args$x[[if (family == "random")
        "method.random.ci" else "method.predict"]]
      for (id in if (shared) 1L else 1:2) {
        parent <- x[[family]][x[[family]]$id == id, ]
        expected <- parent$lower[match(method, parent$method.interval)]
        actual <- if (shared) drawn$args$x[[paste0("lower.", family)]] else
          unname(drawn$args$x[[paste0("lower.", family, ".w")]][id, ])
        expect_equal(actual, expected)
      }
    }
  }
})

test_that("fitted three-level vectors render both shared and component summaries", {
  dat <- data.frame(
    TE = c(-0.4, 0.1, 0.2, 0.7, -0.2, 0.4, 0.8, 1.1, 0.3, 0.9, 1.2, 0.6),
    seTE = rep(c(0.15, 0.2), 6), cluster = rep(1:6, each = 2),
    group = rep(c("A", "B"), each = 6))
  ma <- metagen(TE, seTE, data = dat, cluster = cluster, subgroup = group,
                 sm = "MD", common = FALSE, method.random.ci = "HK",
                 prediction = TRUE, method.predict = "V", tau.common = FALSE)
  expect_length(ma$tau2, 2L)
  for (shared in c(FALSE, TRUE)) {
    x <- metabind(ma, ma, keepmeta = TRUE)
    if (!shared)
      x$data <- list(x$data, x$data)
    before <- x
    expect_silent(capture.output(print(x, overall = TRUE,
      print = list(tau2 = TRUE, tau = TRUE, tau2.ci = TRUE, tau.ci = TRUE))))
    expect_silent(drawn <- interval_slot_draw(x,
      print = list(tau2 = TRUE, tau = TRUE, tau2.ci = TRUE, tau.ci = TRUE)))
    if (shared)
      expect_identical(drawn$args$x$tau2, ma$tau2)
    else
      expect_equal(unname(drawn$args$x$tau2.w), rep(sum(ma$tau2), 2))
    expect_identical(x, before)
  }
})

test_that("custom decreasing transforms retain each interval payload and bound order", {
  fun <- ".meta_interval_slot_decreasing"
  assign(fun, function(x, offset) offset - x, envir = .GlobalEnv)
  on.exit(rm(list = fun, envir = .GlobalEnv))
  x <- interval_slot_fixture(TRUE, FALSE)
  x$overall$func.backtransf <- fun
  x$overall$args.backtransf <- list(list(offset = 3), list(offset = 5))
  before <- x
  expect_silent(output <- capture.output(print(x, overall = TRUE)))
  expect_true(any(grepl("offset = 3", output, fixed = TRUE)))
  expect_silent(drawn <- interval_slot_draw(x))
  native <- drawn$args$x
  for (family in c("random", "predict")) for (id in 1:2) {
    rows <- x[[family]]$id == id
    expect_equal(unname(native[[paste0("lower.", family, ".w")]][id, ]),
                 (1 + 2 * id) - x[[family]]$upper[rows])
    expect_equal(unname(native[[paste0("upper.", family, ".w")]][id, ]),
                 (1 + 2 * id) - x[[family]]$lower[rows])
  }
  expect_identical(x, before)
})
