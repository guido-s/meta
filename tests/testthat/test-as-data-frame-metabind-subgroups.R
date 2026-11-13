metabind_subgroup_fixture <- function() {
  metabin(1:9, 101:109, 2:10, 201:209,
           subgroup = rep(c("A", "B", "C"), each = 3),
           method.random.ci = rep("HK", 3),
           adhoc.hakn.ci = c("", "se", "ci"),
           prediction = TRUE,
           method.predict = c("V", "HK", "HK", "HK-PR"),
           adhoc.hakn.pi = c("", "", "se", "se"))
}

test_that("subgroup random interval keys follow the source matrix columns", {
  ma <- metabind_subgroup_fixture()
  mb <- metabind(ma, ma)
  n.by <- length(ma$subgroup.levels)
  stored <- mb$subgroup$random[mb$subgroup$random$id == 1, ]

  expect_identical(stored$method.interval,
                    rep(unname(ma$method.random.ci), each = n.by))
  expect_identical(stored$adhoc.hakn,
                    rep(unname(ma$adhoc.hakn.ci), each = n.by))
  expect_identical(stored$subgroup.levels,
                    rep(ma$subgroup.levels, length(ma$method.random.ci)))
  expect_equal(as.numeric(stored$lower), as.numeric(ma$lower.random.w))
  expect_equal(stored$Q.b, rep(unname(ma$Q.b.random), each = n.by))

  dat <- as.data.frame(mb, common = FALSE, prediction = FALSE,
                       drop_all_NA = FALSE)
  expect_equal(nrow(dat), 2 * n.by * length(ma$method.random.ci))
  expect_true(all(dat$type == "Random effects model"))
  expect_false(anyNA(dat$TE))
  first <- dat[dat$id == 1, ]
  expect_equal(as.numeric(first$lower.w), as.numeric(ma$lower.random.w))
  expect_equal(first$lower, rep(unname(ma$lower.random), each = n.by))
})

test_that("subgroup prediction keys join to the corresponding pooled intervals", {
  ma <- metabind_subgroup_fixture()
  mb <- metabind(ma, ma)
  n.by <- length(ma$subgroup.levels)
  stored <- mb$subgroup$predict[mb$subgroup$predict$id == 1, ]
  expect_identical(stored$method.interval,
                    rep(unname(ma$method.predict), each = n.by))
  expect_identical(stored$adhoc.hakn,
                    rep(unname(ma$adhoc.hakn.pi), each = n.by))
  expect_equal(as.numeric(stored$lower), as.numeric(ma$lower.predict.w))

  dat <- as.data.frame(mb, common = FALSE, random = FALSE,
                       prediction = TRUE, drop_all_NA = FALSE)
  expect_equal(nrow(dat), 2 * n.by * length(ma$method.predict))
  expect_true(all(dat$type == "Prediction interval"))
  expect_false(anyNA(dat$lower))
  first <- dat[dat$id == 1, ]
  expect_equal(as.numeric(first$lower.w), as.numeric(ma$lower.predict.w))
  expect_equal(first$lower, rep(unname(ma$lower.predict), each = n.by))
})

test_that("subgroup conversion retains joined rows and honors model selection", {
  ma <- metabind_subgroup_fixture()
  mb <- metabind(ma, ma)
  expected <- sum(vapply(mb$subgroup, nrow, integer(1)))
  dat <- as.data.frame(mb)
  expect_equal(nrow(dat), expected)
  expect_false(anyNA(dat$type))
  expect_true("subgroup.levels" %in% names(dat))

  common <- as.data.frame(mb, random = FALSE, prediction = FALSE)
  expect_equal(nrow(common), nrow(mb$subgroup$common))
  expect_true(all(common$type == "Common effect model"))

  overall <- as.data.frame(mb, subgroup = FALSE, common = FALSE,
                           prediction = FALSE)
  expect_equal(nrow(overall), nrow(mb$random))
  expect_true(all(overall$type == "Random effects model"))
  expect_false("subgroup.levels" %in% names(overall))
})

test_that("conversion without source subgroups honors model selection", {
  ma <- metagen(log(2:4), c(0.1, 0.2, 0.3), sm = "RR", prediction = TRUE)
  mb <- metabind(ma, ma)
  random <- as.data.frame(mb, common = FALSE, prediction = FALSE)
  expect_equal(nrow(random), nrow(mb$random))
  expect_true(all(random$type == "Random effects model"))
  prediction <- as.data.frame(mb, common = FALSE, random = FALSE)
  expect_equal(nrow(prediction), nrow(mb$predict))
  expect_true(all(prediction$type == "Prediction interval"))
})

metabind_repeated_interval_fixture <- function(labels = TRUE) {
  ma <- metabin(1:9, 101:109, 2:10, 201:209,
                  subgroup = rep(c("A", "B", "C"), each = 3),
                  method.random.ci = c("HK", "HK"), adhoc.hakn.ci = "se",
                  prediction = TRUE, method.predict = c("V", "V"),
                  adhoc.hakn.pi = "")
  if (labels) {
    ma$text.random <- c("First HK", "Second HK")
    ma$text.predict <- c("First PI", "Second PI")
  }
  metabind(ma, ma, keepmeta = TRUE)
}

test_that("repeated identical methods join by their source interval slots", {
  x <- metabind_repeated_interval_fixture()
  before <- x
  n.by <- length(x$meta[[1]]$subgroup.levels)
  for (model in c("random", "predict")) {
    pooled <- x[[model]]
    subgroup <- x$subgroup[[model]]
    expect_identical(pooled$interval.slot, rep(1:2, 2))
    expect_identical(subgroup$interval.slot, rep(rep(1:2, each = n.by), 2))
    expect_identical(anyDuplicated(pooled[c("id", "interval.slot")]), 0L)
    dat <- as.data.frame(x, common = FALSE, random = model == "random",
                         prediction = model == "predict", drop_all_NA = FALSE)
    expect_equal(nrow(dat), nrow(subgroup))
    expect_false(anyNA(dat$type))
    for (id in 1:2) for (slot in 1:2) {
      rows <- dat$id == id & dat$interval.slot == slot
      parent <- pooled[pooled$id == id & pooled$interval.slot == slot, ]
      expect_equal(sum(rows), n.by)
      expect_identical(dat$text[rows], rep(parent$text, n.by))
      expect_equal(as.numeric(dat$lower[rows]), rep(as.numeric(parent$lower), n.by))
      expect_equal(as.numeric(dat$upper[rows]), rep(as.numeric(parent$upper), n.by))
      limits <- x$meta[[id]][[paste0("lower.", model, ".w")]]
      expect_equal(as.numeric(dat$lower.w[rows]), as.numeric(limits[, slot]))
    }
    if (model == "random")
      expect_identical(dat$df.Q.b, subgroup$df.Q.b)
  }
  expect_identical(x, before)
})

test_that("slot joins survive pooled and subgroup row reordering", {
  x <- metabind_repeated_interval_fixture()
  baseline <- as.data.frame(x, drop_all_NA = FALSE)
  x$random <- x$random[c(4, 2, 1, 3), ]
  x$predict <- x$predict[c(3, 1, 4, 2), ]
  expect_identical(as.data.frame(x, drop_all_NA = FALSE), baseline)
  for (model in c("random", "predict")) {
    x$subgroup[[model]] <- x$subgroup[[model]][rev(seq_len(nrow(x$subgroup[[model]]))), ]
    dat <- as.data.frame(x, common = FALSE, random = model == "random",
                         prediction = model == "predict", drop_all_NA = FALSE)
    for (row in seq_len(nrow(dat))) {
      parent <- x[[model]][x[[model]]$id == dat$id[row] &
                            x[[model]]$interval.slot == dat$interval.slot[row], ]
      expect_identical(dat$text[row], parent$text)
      expect_equal(as.numeric(dat$lower[row]), as.numeric(parent$lower))
    }
  }
})

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
      x$data <- list(x$data, x$data)
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
