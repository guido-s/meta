forest_metabind_subgroup_fixture <- function(common = FALSE, random = TRUE,
                                             keepdata = TRUE, ...) {
  dat <- data.frame(event.e = c(10, 12, 14), n.e = rep(100, 3),
                    event.c = c(10, 12, 14), n.c = rep(100, 3),
                    group1 = c("A", "B", "A"),
                    group2 = c("C", "C", "D"))
  ma <- metabin(event.e, n.e, event.c, n.c, data = dat,
                subgroup = group1, method = "Inverse", method.tau = "DL",
                common = common, random = random, prediction = FALSE, ...)
  metabind(ma, update(ma, subgroup = group2),
           name = c("First grouping", "Second grouping"),
           keepdata = keepdata, keepmeta = TRUE)
}

expected_metabind_data <- function(dat, shared = FALSE) {
  if (shared)
    dat <- dplyr::filter(dat, id == dat$id[1])
  drop <- c("common", "random", "prediction", if (shared) c("id", "name"))
  dat[drop] <- NULL
  dat
}

test_that("shared overall summaries retain all source subgroup IDs", {
  x <- forest_metabind_subgroup_fixture()
  dat <- as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
  subdat <- as.data.frame(x, drop_all_NA = FALSE)
  res <- meta:::prepare_metabind_subgroups(x)

  expect_named(res, c("dat", "subdat", "same.data", "same.results",
                      "shared.overall"))
  expect_true(res$same.data)
  expect_true(res$same.results)
  expect_true(res$shared.overall)
  expect_equal(nrow(dat), 2L)
  expect_equal(nrow(res$dat), 1L)
  expect_equal(res$dat, expected_metabind_data(dat, shared = TRUE))
  expect_false(any(c("id", "name") %in% names(res$dat)))
  expect_equal(res$subdat, subdat)
  expect_equal(nrow(res$subdat), 4L)
  expect_identical(res$subdat$id, subdat$id)
  expect_identical(as.character(res$subdat$subgroup.levels),
                   c("A", "B", "C", "D"))
  expect_equal(meta:::prepare_metabind_subgroups(x, dat = dat), res)
})

test_that("different datasets never share numerically equal summaries", {
  x <- forest_metabind_subgroup_fixture()
  other <- x$data
  other$event.e[1] <- other$event.e[1] + 1
  x$data <- list(x$data, other)
  dat <- as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
  res <- meta:::prepare_metabind_subgroups(x)

  expect_false(res$same.data)
  expect_false(res$same.results)
  expect_false(res$shared.overall)
  expect_equal(res$dat, expected_metabind_data(dat))
  expect_equal(nrow(res$dat), 2L)
  expect_identical(res$dat$id, dat$id)
  expect_equal(res$subdat, as.data.frame(x, drop_all_NA = FALSE))
})

test_that("zero heterogeneity does not collapse different tau methods", {
  x <- forest_metabind_subgroup_fixture()
  x <- metabind(x$meta[[1]], update(x$meta[[2]], method.tau = "REML"))
  dat <- as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
  expect_equal(dat$tau2, c(0, 0))
  expect_equal(dat$TE, c(0, 0))
  expect_equal(dat$lower[1], dat$lower[2])
  expect_equal(dat$upper[1], dat$upper[2])

  res <- meta:::prepare_metabind_subgroups(x)
  expect_true(res$same.data)
  expect_false(res$same.results)
  expect_false(res$shared.overall)
  expect_equal(res$dat, expected_metabind_data(dat))
  expect_identical(res$dat$method.tau, c("DL", "REML"))
})

test_that("common summaries compare original rather than masked tau settings", {
  x <- forest_metabind_subgroup_fixture(common = TRUE, random = FALSE)
  x$overall$method.tau[2] <- "REML"
  dat <- as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
  expect_true(all(is.na(dat$method.tau)))

  res <- meta:::prepare_metabind_subgroups(x)
  expect_true(res$same.data)
  expect_false(res$same.results)
  expect_false(res$shared.overall)
  expect_equal(res$dat, expected_metabind_data(dat))
  expect_equal(nrow(res$dat), 2L)
})

test_that("overall result sets ignore order but retain distinct CI methods", {
  x <- forest_metabind_subgroup_fixture(method.random.ci = c("classic", "HK"))
  x$random <- x$random[c(1, 2, 4, 3), , drop = FALSE]
  dat <- as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
  subdat <- as.data.frame(x, drop_all_NA = FALSE)
  res <- meta:::prepare_metabind_subgroups(x)

  expect_true(res$same.data)
  expect_true(res$same.results)
  expect_true(res$shared.overall)
  expect_equal(nrow(dat), 4L)
  expect_equal(nrow(res$dat), 2L)
  expect_equal(res$dat, expected_metabind_data(dat, shared = TRUE))
  expect_identical(res$dat$method.interval, c("classic", "HK"))
  expect_equal(res$subdat, subdat)
  expect_equal(nrow(res$subdat), 8L)
})

test_that("identical numbers do not erase CI methods within an ID", {
  x <- forest_metabind_subgroup_fixture(method.random.ci = c("classic", "HK"))
  values <- setdiff(names(x$random), c("id", "method.interval"))
  for (id in unique(x$random$id)) {
    rows <- which(x$random$id == id)
    x$random[rows[2], values] <- x$random[rows[1], values]
  }
  res <- meta:::prepare_metabind_subgroups(x)
  expect_true(res$shared.overall)
  expect_equal(nrow(res$dat), 2L)
  expect_identical(res$dat$method.interval, c("classic", "HK"))
  expect_equal(res$dat$TE[1], res$dat$TE[2])
  expect_equal(res$dat$lower[1], res$dat$lower[2])
  expect_equal(nrow(res$subdat), 8L)
})

test_that("selected result sets must match in their entirety", {
  x <- forest_metabind_subgroup_fixture(method.random.ci = c("classic", "HK"))
  x$random <- x$random[-4, , drop = FALSE]
  x$subgroup$random <- subset(x$subgroup$random,
                              !(id == 2 & method.interval == "HK"))
  dat <- as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
  res <- meta:::prepare_metabind_subgroups(x)

  expect_true(res$same.data)
  expect_false(res$same.results)
  expect_false(res$shared.overall)
  expect_equal(res$dat, expected_metabind_data(dat))
  expect_equal(nrow(res$dat), 3L)

  selected <- dat[dat$method.interval == "classic", , drop = FALSE]
  res <- meta:::prepare_metabind_subgroups(x, dat = selected)
  expect_true(res$same.results)
  expect_true(res$shared.overall)
  expect_equal(nrow(res$dat), 1L)
})

test_that("exact duplicate parents never multiply source subgroup rows", {
  for (shared in c(TRUE, FALSE)) {
    x <- forest_metabind_subgroup_fixture()
    if (!shared)
      x$data <- list(x$data, x$data)
    x$subgroup$random <- rbind(x$subgroup$random, x$subgroup$random[1, ])
    expected <- as.data.frame(x, drop_all_NA = FALSE)
    dat <- as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE)
    x$random <- rbind(x$random, x$random[1, ], x$random[2, ])

    res <- meta:::prepare_metabind_subgroups(x)
    expect_identical(res$shared.overall, shared)
    expect_identical(res$same.results, shared)
    expect_equal(nrow(res$dat), if (shared) 1L else 2L)
    expect_equal(res$dat,
                 expected_metabind_data(dat, shared = shared))
    if (!shared)
      expect_identical(res$dat$id, dat$id)
    expect_equal(res$subdat, expected)
    expect_equal(nrow(res$subdat), 5L)
    expect_identical(res$subdat$id, expected$id)
    expect_equal(as.list(res$subdat[1, ]), as.list(res$subdat[5, ]))
  }
})

test_that("all actual metadata and result differences prevent sharing", {
  changes <- list(null.effect = 1, sm = "OR", method.I2 = "tau2",
                   level.ma = 0.99, level.predict = 0.99,
                   level.hetstat = 0.99, tau.preset = 0.1,
                   three.level = TRUE, rho = 0.5, Q = 0.1, tau2 = 0.1,
                   pscale = 100, irscale = 1000, irunit = "years",
                   func.transf = "log", func.backtransf = "exp")
  for (field in names(changes)) {
    x <- forest_metabind_subgroup_fixture()
    if (is.null(x$overall[[field]]))
      x$overall[[field]] <- rep(NA, nrow(x$overall))
    x$overall[[field]][2] <- changes[[field]]
    res <- meta:::prepare_metabind_subgroups(x)
    expect_false(res$same.results, info = field)
    expect_false(res$shared.overall, info = field)
    expect_equal(nrow(res$dat), 2L, info = field)
  }

  x <- forest_metabind_subgroup_fixture()
  x$overall$args.backtransf <- list(list(base = 10), list(base = 2))
  expect_false(meta:::prepare_metabind_subgroups(x)$shared.overall)

  x <- forest_metabind_subgroup_fixture()
  x$random$TE[2] <- x$random$TE[1] + .Machine$double.eps
  res <- meta:::prepare_metabind_subgroups(x)
  expect_false(res$same.results)
  expect_false(res$shared.overall)
  expect_equal(nrow(res$dat), 2L)
})

test_that("raw summary and subgroup columns are neither dropped nor formatted", {
  x <- forest_metabind_subgroup_fixture()
  x$overall$all.missing <- NA_real_
  res <- meta:::prepare_metabind_subgroups(x)
  expect_true("all.missing" %in% names(res$dat))
  expect_true("all.missing" %in% names(res$subdat))
  expect_identical(res$dat$all.missing, NA_real_)
  expect_true(all(is.na(res$subdat$all.missing)))
  expect_type(res$dat$TE, "double")
  expect_type(res$dat$pval, "double")
  expect_type(res$subdat$TE.w, "double")
  expect_equal(res$subdat, as.data.frame(x, drop_all_NA = FALSE))
})

test_that("missing or placeholder datasets cannot enable shared summaries", {
  x <- forest_metabind_subgroup_fixture(keepdata = FALSE)
  expect_null(x$data)
  for (data in list(NULL, data.frame(.not_a_dataset = 1))) {
    x$data <- data
    res <- meta:::prepare_metabind_subgroups(x)
    expect_false(res$same.data)
    expect_false(res$same.results)
    expect_false(res$shared.overall)
    expect_equal(nrow(res$dat), 2L)
    expect_equal(nrow(res$subdat), 4L)
  }
})
