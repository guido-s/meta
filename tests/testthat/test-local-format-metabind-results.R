metabind_format_fixture <- function() {
  data.frame(
    id = rep(1:2, 3), name = rep(c("First", "Second"), 3),
    type = rep(c("Common effect model", "Random effects model",
                 "Prediction interval"), each = 2),
    func = "metagen", sm = "MD", method = "Inverse",
    method.interval = c("classic", "classic", "classic", "HK", "V", "HK"),
    adhoc.hakn = c(NA, NA, NA, "se", NA, "se"),
    level.ma = 0.95, level.predict = 0.90,
    TE = 0.25, seTE = 0.1, lower = 0.05, upper = 0.45,
    statistic = 2.5, pval = 0.0124,
    df = c(NA, NA, Inf, 4, 4, 4))
}

format_metabind_test_results <- function(
    dat, model, print = list(stat = TRUE, pval = TRUE, se = TRUE,
                             func = TRUE, method.interval = TRUE,
                             pval.Q.b = FALSE), ...) {
  meta:::format_metabind_results(
    dat, model, print = print,
    digits = list(digits = 2, se = 3, stat = 2, pval = 4, Q = 2, pval.Q = 4),
    sm.lab = "MD", lab.TE = "TE", ...)
}

test_that("metabind result formatting retains model-specific columns", {
  dat <- metabind_format_fixture()
  common <- format_metabind_test_results(dat, "common")
  random <- format_metabind_test_results(dat, "random")
  prediction <- format_metabind_test_results(dat, "prediction")

  expect_true(all(c("MD", "95% CI", "z", "p-value", "SE(TE)") %in% names(common)))
  expect_false("df" %in% names(common))
  expect_true(all(c("MD", "95% CI", "z|t", "df", "adhoc.hakn") %in% names(random)))
  expect_identical(random$interval, c("classic", "HK"))
  expect_true(all(c("90% PI", "df", "adhoc.hakn") %in% names(prediction)))
  expect_false(any(c("MD", "z", "t", "z|t", "p-value") %in% names(prediction)))
  expect_true(all(common$MD == "0.25"))
  expect_true(all(vapply(list(common, random, prediction), nrow, integer(1)) == 2))
})

test_that("metabind formatting can retain additional columns", {
  dat <- metabind_format_fixture()
  dat$text <- paste("Label", seq_len(nrow(dat)))
  for (model in c("common", "random", "prediction")) {
    default <- format_metabind_test_results(dat, model)
    retained <- format_metabind_test_results(dat, model, keep = c("type", "text"))
    rows <- match(retained$text, dat$text)

    expect_false(any(c("type", "text") %in% names(default)))
    expect_identical(retained$type, dat$type[rows])
    expect_identical(retained$text, dat$text[rows])
    expect_identical(retained[, names(default), drop = FALSE], default)
  }
})

test_that("normalized pooled data use the same metabind formatter", {
  dat <- metabind_format_fixture()
  pooled <- dat
  fields <- c("TE", "seTE", "lower", "upper", "statistic", "pval", "df")
  for (field in fields)
    names(pooled)[names(pooled) == field] <- paste0(field, ".pooled")
  pooled$id <- pooled$name <- NULL
  names(pooled) <- sub("\\.pooled$", "", names(pooled))
  dat$id <- dat$name <- NULL

  for (model in c("common", "random", "prediction"))
    expect_identical(format_metabind_test_results(pooled, model),
                     format_metabind_test_results(dat, model))
})

test_that("metabind formatting accepts shared interval headings and level visibility", {
  dat <- metabind_format_fixture()
  tables <- lapply(c("common", "random", "prediction"), function(model)
    format_metabind_test_results(dat, model, lab.interval = "CI/PI",
                                 show.level = TRUE))
  combined <- dplyr::bind_rows(tables)
  expect_true("CI/PI" %in% names(combined))
  expect_false(any(c("95% CI", "90% PI") %in% names(combined)))
  expect_equal(trimws(combined$level), c(rep("0.9500", 4), rep("0.9000", 2)))

  dat$level.predict <- 0.95
  for (model in c("common", "random", "prediction")) {
    formatted <- format_metabind_test_results(
      dat, model, lab.interval = "95% CI/PI", show.level = FALSE)
    expect_true("95% CI/PI" %in% names(formatted))
    expect_false("level" %in% names(formatted))
  }
})

test_that("metabind formatting drops all-NA columns after column selection", {
  dat <- metabind_format_fixture()
  dat$func <- NA_character_
  dat$statistic <- dat$pval <- dat$df <- NA_real_
  dat$subgroup.levels <- rep(c("A", "B"), 3)
  dat$Q.b <- NA_real_
  formatted <- format_metabind_test_results(dat, "random", with.subgroups = TRUE)

  expect_false(any(c("func", "statistic", "z|t", "p-value", "df", "Q.b") %in%
                     names(formatted)))
  expect_identical(formatted$subgroup.levels, c("A", "B"))
  expect_false(any(vapply(formatted, function(column) all(is.na(column)), logical(1))))
})

test_that("subgroup tests retain one result per component and method before hiding methods", {
  dat <- metabind_format_fixture()[rep(c(3, 4, 3, 4), each = 2), ]
  dat$id <- rep(c(1, 1, 2, 2), each = 2)
  dat$Q.b <- rep(c(0.375, 0.0778, 0.2, 0.1), each = 2)
  dat$df.Q.b <- rep(list(1), nrow(dat))
  dat$test.Q.b <- "Q"
  dat$pval.Q.b <- rep(c(0.5403, 0.7803, 0.6547, 0.7518), each = 2)
  flags <- list(stat = TRUE, pval = TRUE, se = FALSE, method.interval = FALSE,
                 pval.Q.b = TRUE)
  for (ad.hoc in c(FALSE, TRUE)) {
    if (ad.hoc) {
      dat$method.interval <- "HK"
      dat$adhoc.hakn <- rep(c("", "se", "", "se"), each = 2)
    }
    formatted <- format_metabind_test_results(
      dat, "random", print = flags, with.subgroups = TRUE)
    first <- c(1, 3, 5, 7)
    expect_equal(trimws(formatted[["Q.b"]][first]),
                  trimws(meta:::formatN(dat$Q.b[first], digits = 2)))
    expect_true(all(trimws(formatted[["Q.b"]][-first]) == "."))
    expect_true(all(trimws(formatted[["d.f."]][first]) == "1"))
    expect_true(all(trimws(formatted[["d.f."]][-first]) == "."))
    expect_equal(trimws(formatted[["p-value"]][first]),
                  trimws(meta:::formatPT(dat$pval.Q.b[first], digits = 4)))
    expect_true(all(formatted$test.Q.b[first] == "Q"))
    expect_true(all(is.na(formatted$test.Q.b[-first])))
    expect_false(any(c("interval", "adhoc.hakn", "df", "df.Q.b", "pval.Q.b") %in%
                       names(formatted)))
    expect_false(anyDuplicated(names(formatted)) > 0)
  }
})

test_that("identical method occurrences retain their own subgroup test rows", {
  dat <- metabind_format_fixture()[rep(4, 4), ]
  dat$interval.slot <- rep(1:2, each = 2)
  dat$subgroup.levels <- rep(c("A", "B"), 2)
  dat$Q.b <- rep(c(0.25, 0.75), each = 2)
  dat$df.Q.b <- rep(list(c(1, 6.42)), 4)
  dat$test.Q.b <- "F"
  dat$pval.Q.b <- rep(c(0.6, 0.4), each = 2)
  formatted <- format_metabind_test_results(
    dat, "random", with.subgroups = TRUE,
    print = list(stat = TRUE, pval = TRUE, se = FALSE, pval.Q.b = TRUE))
  expect_identical(trimws(formatted$Q.b), c("0.25", ".", "0.75", "."))
  expect_identical(trimws(formatted[["d.f."]]), c("1, 6.42", ".", "1, 6.42", "."))
  expect_false("interval.slot" %in% names(formatted))
})

test_that("unavailable subgroup tests remain missing without requiring raw columns", {
  dat <- metabind_format_fixture()
  flags <- list(stat = TRUE, pval = TRUE, se = FALSE, pval.Q.b = TRUE)
  for (model in c("common", "random")) {
    formatted <- format_metabind_test_results(
      dat, model, print = flags, with.subgroups = TRUE)
    expect_true(all(trimws(formatted[["Q.b"]]) == "."))
    expect_true(all(trimws(formatted[["d.f."]]) == "."))
    expect_true(all(trimws(formatted[["p-value"]]) == "."))
  }
})

test_that("numeric result columns honor missing labels and thousands separators", {
  dat <- metabind_format_fixture()
  dat$k <- rep(c(1000, NA_real_), 3)
  dat$df[3:6] <- rep(c(1200.126, NA_real_), 2)
  dat$level.ma <- dat$level.predict <- rep(c(0.95, NA_real_), 3)
  for (model in c("common", "random", "prediction")) {
    formatted <- format_metabind_test_results(
      dat, model, show.level = TRUE, lab.NA = "MISSING", big.mark = ",")
    expect_identical(trimws(formatted$k), c("1,000", "MISSING"))
    expect_identical(trimws(formatted$level), c("0.9500", "MISSING"))
    if (model != "common")
      expect_identical(trimws(formatted$df), c("1,200.13", "MISSING"))
  }

  formatted <- format_metabind_test_results(
    dat, "random", with.subgroups = TRUE, show.level = TRUE,
    lab.NA = "MISSING", big.mark = ",")
  expect_identical(formatted$level, c(0.95, NA_real_))
  expect_identical(trimws(formatted$k), c("1,000", "MISSING"))
})

test_that("entirely missing numeric result columns are still removed", {
  dat <- metabind_format_fixture()
  dat$k <- dat$df <- dat$level.ma <- dat$level.predict <- NA_real_
  for (model in c("common", "random", "prediction")) {
    formatted <- format_metabind_test_results(
      dat, model, show.level = TRUE, lab.NA = "MISSING")
    expect_false(any(c("k", "df", "level") %in% names(formatted)))
  }
})
