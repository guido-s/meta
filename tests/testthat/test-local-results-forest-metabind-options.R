test_that("mixed rate transformations use time harmonic means and rate scaling", {
  ma <- metarate(c(2, 4, 8), c(100, 200, 300), sm = "IRLN", irscale = 100)
  other <- update(ma, sm = "IRFT")
  res <- forest_metabind_options(metabind(ma, other))
  expect_equal(res$data$TE[res$data$id == 1],
               rep(meta:::backtransf(ma$TE.common, ma$sm) * 100, 2))
  expect_equal(res$data$TE[res$data$id == 2],
               rep(meta:::backtransf(other$TE.common, other$sm,
                                    time = other$t.harmonic.mean) * 100, 2))
})

test_that("prediction rows retain tau estimates and methods from selected random results", {
  ma <- metagen_forest_options(prediction = TRUE, method.predict = c("V", "HK"))
  other <- update(ma, method.tau = "REML", random = FALSE)
  mb <- metabind(ma, other)
  res <- forest_metabind_options(
    mb, lab.NA = ".", print = list(tau2 = TRUE, tau = TRUE,
                                   tau2.ci = TRUE, tau.ci = TRUE,
                                   method.tau = TRUE))
  dat <- res$data
  random.row <- which(dat$id == 1 & dat$type == "Random effects model")[1]
  prediction.rows <- which(dat$id == 1 & dat$type == "Prediction interval")
  for (nm in c("tau2", "tau", "tau2.ci", "tau.ci", "method.tau"))
    expect_equal(dat[[nm]][prediction.rows],
                 rep(dat[[nm]][random.row], length(prediction.rows)))
  expect_true(all(dat$method.tau[prediction.rows] == "DL"))

  without.random <- dat$id == 2 & dat$type == "Prediction interval"
  expect_true(all(dat$tau2[without.random] == "."))
  expect_true(all(dat$method.tau[without.random] == "."))
  hidden <- forest_metabind_options(mb, random = FALSE, lab.NA = ".",
                                    print = list(tau2 = TRUE, tau = TRUE,
                                                 method.tau = TRUE))$data
  expect_true(all(hidden$tau2[hidden$type == "Prediction interval"] == "."))
  expect_true(all(hidden$method.tau[hidden$type == "Prediction interval"] == "."))
})

test_that("custom transformation metadata keeps named arguments and missing functions", {
  fun.name <- ".meta_forest_test_shift"
  assign(fun.name, function(x, offset = 0, scale = 1) (x + offset) * scale,
         envir = .GlobalEnv)
  on.exit(rm(list = fun.name, envir = .GlobalEnv))
  ma <- metagen(1:3, rep(1, 3), sm = "MD", func.backtransf = fun.name,
                 args.backtransf = list(scale = 100))
  other <- metagen(1:3, rep(1, 3), sm = "MD")
  res <- forest_metabind_options(metabind(ma, other))
  expect_equal(res$data$TE[res$data$id == 1],
               res$data$TE.orig[res$data$id == 1] * 100)
  expect_equal(res$data$TE[res$data$id == 2],
               res$data$TE.orig[res$data$id == 2])
})

test_that("relative measures respect non-exponential custom transformations", {
  ma <- metagen(c(0.2, 0.3, 0.4), rep(0.1, 3), sm = "RR",
                 func.backtransf = "plogis", prediction = TRUE)
  res <- forest_metabind_options(metabind(ma, ma), digits = 4)
  expect_equal(res$data$TE, plogis(res$data$TE.orig))
  expect_equal(res$data$lower, plogis(res$data$lower.orig))
  expect_equal(res$data$upper, plogis(res$data$upper.orig))
  expect_true("0.5744" %in% trimws(res$effect.format))

  unchanged <- forest_metabind_options(metabind(ma, ma), backtransf = FALSE)
  expect_equal(unchanged$data$TE, unchanged$data$TE.orig)
  expect_equal(unchanged$data$lower, unchanged$data$lower.orig)
  expect_equal(unchanged$data$upper, unchanged$data$upper.orig)
})

test_that("mixed effective transformations disable the shared exponential shortcut", {
  custom <- metagen(c(0.2, 0.3, 0.4), rep(0.1, 3), sm = "RR",
                     func.backtransf = "plogis")
  standard <- metagen(c(0.2, 0.3, 0.4), rep(0.1, 3), sm = "OR")
  mixed <- forest_metabind_options(metabind(custom, standard))$data
  expect_equal(mixed$TE[mixed$id == 1], plogis(mixed$TE.orig[mixed$id == 1]))
  expect_equal(mixed$TE[mixed$id == 2], exp(mixed$TE.orig[mixed$id == 2]))

  standard <- forest_metabind_options(metabind(standard, standard), digits = 4)
  expect_equal(standard$data$TE, standard$data$TE.orig)
  expect_true("1.3499" %in% trimws(standard$effect.format))
})
