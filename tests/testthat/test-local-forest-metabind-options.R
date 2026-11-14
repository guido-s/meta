test_that("constant metadata columns are automatic and print flags override", {
  ma <- metagen_forest_options(common = TRUE, random = FALSE,
                              prediction = FALSE)
  mb <- metabind(ma, ma)
  cols <- c("func", "method", "method.interval", "sm")
  flags <- setNames(rep(list(TRUE), length(cols)), cols)

  res <- forest_metabind_options(mb)
  expect_s3_class(res$data, "data.frame")
  expect_true(all(c("id", "type", cols) %in% names(res$data)))
  expect_equal(length(res$leftlabs), length(res$leftcols))
  expect_false(any(cols %in% res$leftcols))
  expect_true(all(cols %in% forest_metabind_options(mb, print = flags)$leftcols))

  flags[] <- list(FALSE)
  expect_false(any(cols %in%
                   forest_metabind_options(mb, print = flags)$leftcols))
})

test_that("print.method.tau controls the heterogeneity estimator column", {
  ma <- metagen_forest_options()
  mb <- metabind(ma, update(ma, method.tau = "REML"))
  expect_true("method.tau" %in% forest_metabind_options(mb)$leftcols)
  expect_false("method.tau" %in%
                 forest_metabind_options(mb, print.method.tau = FALSE)$leftcols)
  expect_false("method.tau" %in%
                 forest_metabind_options(mb, print = list(method.tau = FALSE))$leftcols)
  expect_true("method.tau" %in%
                forest_metabind_options(metabind(ma, ma),
                                        print.method.tau = TRUE)$leftcols)
  expect_true("method.tau" %in%
                forest_metabind_options(mb, print = list(method = FALSE,
                                                        method.tau = TRUE))$leftcols)
})

test_that("automatic I2 columns follow tau estimator methods", {
  data("Fleiss1993cont", package = "meta")
  ma <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
                  data = Fleiss1993cont, sm = "SMD", common = FALSE,
                  method.tau = "REML")
  mb <- metabind(ma, update(ma, method.tau = "DL"))
  drawn <- forest_metabind_options(
    mb, print = list(I2 = TRUE, method.I2 = TRUE, I2.ci = FALSE))
  expect_true(all(c("method.tau", "I2") %in% drawn$leftcols))
  expect_gt(match("I2", drawn$leftcols), match("method.tau", drawn$leftcols))
  expect_equal(match("method.I2", drawn$leftcols), match("I2", drawn$leftcols) + 1)
  intervals <- forest_metabind_options(
    mb, print = list(I2 = TRUE, I2.ci = TRUE, tau2 = TRUE, tau2.ci = TRUE,
                     tau = TRUE, tau.ci = TRUE, H = TRUE, H.ci = TRUE,
                     Q = TRUE, pval.Q = TRUE, method.I2 = TRUE,
                     method.tau = TRUE, method.tau.ci = TRUE))
  expect_gt(match("I2", intervals$leftcols), match("method.tau.ci", intervals$leftcols))
  expect_equal(match("I2.ci", intervals$leftcols), match("I2", intervals$leftcols) + 1)
  ordered <- c("tau2", "tau2.ci", "tau", "tau.ci", "method.tau",
               "method.tau.ci", "I2", "I2.ci", "method.I2", "H", "H.ci",
               "Q", "df.Q", "pval.Q")
  expect_identical(intervals$leftcols[intervals$leftcols %in% ordered], ordered)
  columns <- c("studlab", "method.I2", "I2", "method.tau", "effect", "ci")
  explicit <- forest_metabind_options(mb, leftcols = columns)
  expect_identical(explicit$leftcols, columns)
})

test_that("automatic metadata columns use every user-defined subgroup", {
  binary <- metabin(c(12, 20, 30), c(100, 100, 100),
                    c(8, 10, 12), c(100, 100, 100), sm = "RR",
                    method = "MH", common = TRUE, random = FALSE)
  generic <- metagen_forest_options(common = FALSE, random = TRUE,
                                   prediction = FALSE,
                                   method.random.ci = "HK")
  mb <- metabind(binary, generic, prediction = FALSE)
  cols <- c("func", "method", "method.interval", "sm")

  ungrouped <- forest_metabind_options(mb)
  expect_true(all(c("func", "method", "method.interval") %in%
                    ungrouped$leftcols))
  expect_false("sm" %in% ungrouped$leftcols)
  grouped <- forest_metabind_options(mb, subgroup = type)
  expect_false(any(cols %in% grouped$leftcols))
  expect_identical(grouped$data$type, ungrouped$data$type)
  expect_identical(grouped$data$id, ungrouped$data$id)
  expect_equal(as.character(grouped$data$type),
               c("Common effect model", "Random effects model"))

  flags <- setNames(rep(list(TRUE), length(cols)), cols)
  expect_true(all(cols %in%
                  forest_metabind_options(mb, subgroup = type,
                                          print = flags)$leftcols))
  flags[] <- list(FALSE)
  expect_false(any(cols %in%
                   forest_metabind_options(mb, print = flags)$leftcols))

  varied <- metabind(binary, generic,
                    update(generic, method.random.ci = "classic"),
                    prediction = FALSE)
  res <- forest_metabind_options(varied, subgroup = type)
  expect_true("method.interval" %in% res$leftcols)
  expect_false(any(c("func", "method", "sm") %in% res$leftcols))
})

test_that("p-value and statistic defaults depend on any nonmissing p-values", {
  ma <- metagen_forest_options(common = TRUE, random = TRUE,
                              prediction = TRUE, method.predict = "V")
  mb <- metabind(ma, ma)
  res <- forest_metabind_options(mb)
  expect_true(any(!is.na(res$data$pval)))
  expect_true(any(is.na(res$data$pval)))
  expect_true(all(c("pval", "statistic") %in% res$leftcols))
  expect_false(any(c("pval", "statistic") %in%
                   forest_metabind_options(
                     mb, print = list(pval = FALSE, stat = FALSE))$leftcols))

  mb$common$pval[] <- NA_real_
  mb$random$pval[] <- NA_real_
  expect_false(any(c("pval", "statistic") %in%
                   forest_metabind_options(mb)$leftcols))
  expect_false(any(c("pval", "statistic") %in%
                   forest_metabind_options(mb, common = FALSE,
                                           random = FALSE)$leftcols))
})

test_that("constant CI and PI levels are included in the interval heading", {
  ma <- metagen_forest_options(common = TRUE, random = TRUE,
                              prediction = TRUE, method.predict = "V",
                              level.ma = 0.95, level.predict = 0.95)
  same <- forest_metabind_options(metabind(ma, ma))
  expect_equal(same$leftlabs[match("ci", same$leftcols)], "95% CI/PI")
  expect_false("level" %in% same$leftcols)

  ma <- update(ma, level.predict = 0.99)
  different <- forest_metabind_options(metabind(ma, ma))
  expect_equal(different$leftlabs[match("ci", different$leftcols)],
               "95% CI/99% PI")
  expect_false("level" %in% different$leftcols)
})

test_that("varying levels of either interval kind require a level column", {
  ma <- metagen_forest_options(common = TRUE, random = TRUE,
                              prediction = TRUE, method.predict = "V",
                              level.ma = 0.95, level.predict = 0.95)
  for (other in list(update(ma, level.ma = 0.90),
                     update(ma, level.predict = 0.99))) {
    res <- forest_metabind_options(metabind(ma, other))
    expect_true("level" %in% res$leftcols)
    expect_true("level" %in% names(res$data))
    expect_equal(res$leftlabs[match("ci", res$leftcols)], "CI/PI")
    expect_equal(length(unique(res$data$level)), 2L)
  }
})

test_that("heterogeneity uses preferred rows and tau is retained for predictions", {
  ma <- metagen_forest_options(common = TRUE, random = TRUE,
                              prediction = TRUE,
                              method.random.ci = c("classic", "HK"),
                              method.predict = c("V", "HK"))
  mb <- metabind(ma, ma)
  cols <- c("tau2", "tau", "I2", "Q", "pval.Q")
  flags <- setNames(rep(list(TRUE), length(cols)), cols)
  selections <- list(
    list(common = TRUE, random = TRUE, prediction = TRUE),
    list(common = TRUE, random = FALSE, prediction = TRUE),
    list(common = FALSE, random = FALSE, prediction = TRUE)
  )
  preferred <- c("Random effects model", "Common effect model",
                 "Prediction interval")

  for (i in seq_along(selections)) {
    res <- do.call(forest_metabind_options,
                   c(list(x = mb, print = flags, lab.NA = ""), selections[[i]]))
    expect_true(all(cols %in% res$leftcols))
    expect_true(all(cols %in% names(res$data)))
    for (id in unique(res$data$id)) {
      rows <- which(res$data$id == id)
      first <- rows[which(res$data$type[rows] == preferred[i])[1]]
      expect_false(is.na(first))
      for (col in cols) {
        values <- res$data[[col]][rows]
        present <- !is.na(values) & trimws(as.character(values)) != ""
        expected <- if (i == 2 && col %in% c("tau2", "tau"))
          integer()
        else if (i == 1 && col %in% c("tau2", "tau"))
          c(first, rows[res$data$type[rows] == "Prediction interval"])
        else
          first
        expect_equal(rows[present], expected,
                      info = paste(preferred[i], "id", id, col))
      }
    }
  }
})

test_that("common effect rows never show between-study variance estimates", {
  first <- metagen(1:3, 3:1 / 10, random = FALSE)
  fields <- c("tau2", "tau", "tau2.ci", "tau.ci", "method.tau", "method.tau.ci")
  for (other.random in c(FALSE, TRUE)) {
    second <- metagen(1:3, 3:1 / 10, random = other.random)
    mb <- metabind(first, second)
    automatic <- forest_metabind_options(mb)
    common <- automatic$data$type == "Common effect model"
    expect_true(all(automatic$data$tau2[common] == meta::gs("lab.NA")))
    expect_true(all(automatic$data$tau[common] == meta::gs("lab.NA")))
    drawn <- forest_metabind_options(
      mb, lab.NA = "MISSING",
      print = list(tau2 = TRUE, tau = TRUE, tau2.ci = TRUE, tau.ci = TRUE,
                   method.tau = TRUE, method.tau.ci = TRUE))
    common <- drawn$data$type == "Common effect model"
    for (field in fields)
      expect_true(all(drawn$data[[field]][common] == "MISSING"), info = field)
    if (other.random) {
      random <- drawn$data$type == "Random effects model"
      expect_true(all(drawn$data$tau2[random] != "MISSING"))
      expect_true(all(drawn$data$tau[random] != "MISSING"))
    }
  }
})

test_that("missing and suppressed heterogeneity values use lab.NA", {
  ma <- metagen_forest_options()
  mb <- metabind(ma, ma)
  mb$overall$tau2[2] <- NA_real_
  mb$overall$method.tau[2] <- NA_character_

  for (grouped in c(FALSE, TRUE)) {
    res <- if (grouped)
      forest_metabind_options(mb, subgroup = type, lab.NA = "MISSING",
                              print = list(method = TRUE, tau2 = TRUE))
    else
      forest_metabind_options(mb, lab.NA = "MISSING",
                              print = list(method = TRUE, tau2 = TRUE))
    selected <- res$data$id == 2 & res$data$type == "Random effects model"
    duplicate <- res$data$id == 2 & res$data$type == "Common effect model"
    expect_equal(res$data$tau2[selected], "MISSING")
    expect_equal(res$data$method.tau[selected], "MISSING")
    expect_equal(res$data$tau2[duplicate], "MISSING")
    expect_equal(res$data$method.tau[duplicate], "MISSING")
  }
})

test_that("empty tau estimator names use lab.NA without changing stored methods", {
  mbin1 <- metabin(1:5, 101:105, 6:10, 201:205, sm = "OR", method = "Inverse")
  mbin2 <- update(mbin1, method = "LRP")
  mb <- metabind(mbin1, mbin2)
  original <- mb$overall$method.tau
  expect_identical(original[2], "")
  for (missing in c("", ".", "MISSING")) {
    res <- forest_metabind_options(mb, lab.NA = missing)
    expect_true("method.tau" %in% res$leftcols)
    lrp <- res$data$id == 2 & res$data$type == "Random effects model"
    inverse <- res$data$id == 1 & res$data$type == "Random effects model"
    expect_true(any(lrp))
    expect_true(all(res$data$method.tau[lrp] == missing))
    expect_true(all(res$data$method.tau[inverse] == mbin1$method.tau))
  }
  expect_identical(mb$overall$method.tau, original)
})

test_that("display skeleton does not pool component result rows", {
  ma <- metagen_forest_options(prediction = TRUE)
  mb <- metabind(ma, ma)
  expect_warning(res <- forest_metabind_options(mb, hetstat = TRUE),
                 "Input to argument 'hetstat' ignored", fixed = TRUE)
  expect_equal(nrow(res$data), nrow(as.data.frame(mb)))
})
