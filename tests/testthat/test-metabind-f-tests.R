metabind_f_test_fixture <- function(methods = "HK", common = FALSE, n.by = 2L) {
  dat <- data.frame(
    event = rep(seq(10, 50, 10), n.by) + rep(20 * (seq_len(n.by) - 1), each = 5),
    total = 100, group = rep(LETTERS[seq_len(n.by)], each = 5))
  metaprop(event, total, data = dat, method = "GLMM", method.random.ci = methods,
             common = common, subgroup = group, tau.common = TRUE)
}

metabind_f_test_forest <- function(x, ...) {
  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off())
  meta::forest(x, layout = "meta", details = FALSE, ...)
}

test_that("subgroup DF formatting preserves pairs, fractional DF, infinity, and missing values", {
  expect_identical(meta:::format_df(c(1, 8)), "1, 8")
  expect_identical(meta:::format_df(c(1, 1)), "1, 1")
  expect_identical(meta:::format_df(c(1, 6.42)), "1, 6.42")
  expect_identical(meta:::format_df(c(1, Inf)), "1, Inf")
  expect_identical(meta:::format_df(c(1000, 2000), big.mark = ","),
                   "1,000, 2,000")
  expect_identical(meta:::format_df(NULL, lab.NA = "MISSING"), "MISSING")
  expect_identical(meta:::format_df(NA_real_, lab.NA = "MISSING"), "MISSING")
  expect_identical(meta:::subgroup_test_heading(c("Q", NA_character_)), "Q")
  expect_identical(meta:::subgroup_test_heading(c("F", "F")), "F")
  expect_identical(meta:::subgroup_test_heading(c("Q", "F")), "Q/F")
  expect_identical(meta:::subgroup_test_heading(NA_character_), "Q/F")
})

test_that("bound subgroup tests recover missing P and preserve stored missing values", {
  ma <- metabind_f_test_fixture(c("HK", "classic"), common = TRUE)
  baseline <- metabind(ma, ma)
  for (model in c("common", "random")) {
    field <- paste0("pval.Q.b.", model)
    missing <- ma
    missing[[field]] <- NULL
    actual <- metabind(missing, missing)
    expect_equal(actual$subgroup[[model]]$pval.Q.b,
                 baseline$subgroup[[model]]$pval.Q.b)
    expect_null(missing[[field]])
    for (value in c(0.4321, NA_real_)) {
      stored <- ma
      stored[[field]] <- rep(value, length(ma[[paste0("Q.b.", model)]]))
      actual <- metabind(stored, stored)
      expect_identical(actual$subgroup[[model]]$pval.Q.b,
                       rep(value, nrow(actual$subgroup[[model]])))
    }
  }
})

test_that("unavailable subgroup tests do not get recycled to later interval slots", {
  ma <- metabind_f_test_fixture(c("HK", "classic"))
  partial <- ma
  partial$Q.b.random <- ma$Q.b.random[1]
  partial$df.Q.b.random <- ma$df.Q.b.random[[1]]
  partial$pval.Q.b.random <- NULL
  bound <- metabind(partial, partial)
  for (id in 1:2) {
    dat <- bound$subgroup$random[bound$subgroup$random$id == id, ]
    first <- dat$interval.slot == 1L
    expect_equal(dat$Q.b[first], rep(unname(ma$Q.b.random[1]), 2))
    expect_equal(dat$pval.Q.b[first], rep(unname(ma$pval.Q.b.random[1]), 2))
    expect_identical(dat$df.Q.b[first], rep(list(unname(ma$df.Q.b.random[[1]])), 2))
    expect_identical(dat$test.Q.b[first], rep("F", 2))
    expect_true(all(is.na(dat$Q.b[!first])))
    expect_true(all(is.na(dat$pval.Q.b[!first])))
    expect_true(all(vapply(dat$df.Q.b[!first], function(df) all(is.na(df)), logical(1))))
    expect_true(all(is.na(dat$test.Q.b[!first])))
  }
  partial$Q.b.random <- numeric()
  partial$df.Q.b.random <- list()
  bound <- metabind(partial, partial)
  expect_equal(nrow(bound$subgroup$random), 8L)
  expect_true(all(is.na(bound$subgroup$random$Q.b)))
  expect_true(all(is.na(bound$subgroup$random$pval.Q.b)))
  expect_true(all(is.na(bound$subgroup$random$test.Q.b)))
})

test_that("bound Q tests retain the generic DF fallback", {
  ma <- metagen(c(0.1, 0.4, 0.2, 0.8, 1.2, 1.5, 0.9, 1.8), rep(0.2, 8),
                 subgroup = rep(c("A", "B"), each = 4), sm = "MD",
                 method.random.ci = c("HK", "classic"))
  baseline <- metabind(ma, ma)
  ma$df.Q.b.common <- ma$df.Q.b.random <- NULL
  ma$pval.Q.b.common <- ma$pval.Q.b.random <- NULL
  actual <- metabind(ma, ma)
  for (model in c("common", "random")) {
    expect_identical(actual$subgroup[[model]]$df.Q.b, baseline$subgroup[[model]]$df.Q.b)
    expect_equal(actual$subgroup[[model]]$pval.Q.b, baseline$subgroup[[model]]$pval.Q.b)
  }
})

test_that("binding retains complete F-test DF without duplicating subgroup rows", {
  for (n.by in c(2L, 3L)) {
    ma <- metabind_f_test_fixture(common = TRUE, n.by = n.by)
    mb <- metabind(ma, ma, name = c("First", "Second"))
    expect_equal(nrow(mb$subgroup$common), 2L * n.by)
    expect_equal(nrow(mb$subgroup$random), 2L * n.by)
    expect_identical(mb$subgroup$random$df.Q.b,
                     rep(list(unname(ma$df.Q.b.random)), 2L * n.by))
    expect_true(all(mb$subgroup$common$test.Q.b == "Q"))
    expect_true(all(mb$subgroup$random$test.Q.b == "F"))
    expect_equal(unname(ma$df.Q.b.random), c(n.by - 1L, 4L * n.by))
    expect_equal(unname(ma$pval.Q.b.random),
                  pf(ma$Q.b.random, ma$df.Q.b.random[1], ma$df.Q.b.random[2],
                     lower.tail = FALSE))

    dat <- as.data.frame(mb, drop_all_NA = FALSE)
    expect_equal(nrow(dat), 4L * n.by)
    expect_type(dat$df.Q.b, "list")
    expect_true(all(lengths(dat$df.Q.b[dat$type == "Random effects model"]) == 2L))
    prepared <- meta:::prepare_metabind_subgroups(mb)
    expect_equal(nrow(prepared$subdat), 4L * n.by)
    expect_identical(prepared$subdat$df.Q.b, dat$df.Q.b)
  }
})

test_that("bound classic and HK tests retain their own metadata in both method orders", {
  for (methods in list(c("classic", "HK"), c("HK", "classic"))) {
    ma <- metabind_f_test_fixture(methods)
    mb <- metabind(ma, ma, name = c("First", "Second"))
    expect_equal(nrow(mb$subgroup$random), 8L)
    for (id in 1:2) {
      dat <- mb$subgroup$random[mb$subgroup$random$id == id, ]
      expect_identical(dat$method.interval, rep(methods, each = 2))
      expect_identical(dat$df.Q.b, rep(unname(ma$df.Q.b.random), each = 2))
      expect_identical(dat$test.Q.b,
                       rep(ifelse(methods == "HK", "F", "Q"), each = 2))
      expect_equal(unname(dat$pval.Q.b), rep(unname(ma$pval.Q.b.random), each = 2))
    }
    output <- capture.output(print(mb, overall = FALSE, details.methods = FALSE))
    expect_equal(sum(grepl("Q/F +d[.]f[.]", output)), 2)
    expect_equal(sum(grepl("1, 8", output, fixed = TRUE)), 2)
    forest <- metabind_f_test_forest(mb)
    expect_identical(forest$leftlabs[match("Q.b", forest$leftcols)], "Q/F")
    first <- forest$data$subgroup.levels == "A"
    expect_identical(forest$data$df.Q.b[first],
                     rep(ifelse(methods == "HK", "1, 8", "1"), 2))
  }
})

test_that("printed headings reflect each component while forest headings reflect all displayed tests", {
  ma <- metabind_f_test_fixture(common = TRUE)
  mb <- metabind(ma, ma, name = c("First", "Second"))
  output <- capture.output(print(mb, overall = FALSE, details.methods = FALSE))
  expect_equal(sum(grepl(" Q +d[.]f[.]", output)), 2)
  expect_equal(sum(grepl(" F +d[.]f[.]", output)), 2)
  expect_equal(sum(grepl("1, 8", output, fixed = TRUE)), 2)
  expect_false(any(grepl("test.Q.b", output, fixed = TRUE)))

  mixed <- metabind_f_test_forest(mb)
  expect_identical(mixed$leftlabs[match("Q.b", mixed$leftcols)], "Q/F")
  random <- metabind_f_test_forest(mb, common = FALSE)
  expect_identical(random$leftlabs[match("Q.b", random$leftcols)], "F")
  common <- metabind_f_test_forest(mb, random = FALSE)
  expect_identical(common$leftlabs[match("Q.b", common$leftcols)], "Q")

  mb <- metabind(update(ma, common = FALSE, method.random.ci = "classic"),
                  update(ma, common = FALSE), name = c("Classic", "HK"))
  output <- capture.output(print(mb, overall = FALSE, details.methods = FALSE))
  expect_equal(sum(grepl(" Q +d[.]f[.]", output)), 1)
  expect_equal(sum(grepl(" F +d[.]f[.]", output)), 1)
  expect_false(any(grepl("Q/F", output, fixed = TRUE)))
})

test_that("bound F-test columns respect overrides and preserve repeated DF values", {
  ma <- metabind_f_test_fixture()
  mb <- metabind(ma, ma)
  explicit <- metabind_f_test_forest(
    mb, leftcols = c("studlab", "test.Q.b", "Q.b", "df.Q.b", "pval.Q.b"))
  expect_identical(explicit$leftlabs[match("test.Q.b", explicit$leftcols)], "Test")
  expect_identical(explicit$data$test.Q.b, c("F", "", "F", ""))
  hidden <- metabind_f_test_forest(mb, print = list(Q.b = FALSE, pval.Q.b = FALSE))
  expect_false(any(c("Q.b", "df.Q.b", "pval.Q.b", "test.Q.b") %in% hidden$leftcols))
  output <- capture.output(print(mb, overall = FALSE,
                                  print = list(pval.Q.b = FALSE), details.methods = FALSE))
  expect_false(any(grepl("d[.]f[.]|Q[.]b|test[.]Q[.]b", output)))

  mb$subgroup$random$df.Q.b <- rep(list(c(1, 1)), 4)
  forest <- metabind_f_test_forest(mb)
  expect_identical(forest$data$df.Q.b, c("1, 1", "", "1, 1", ""))
  mb$subgroup$random$df.Q.b <- rep(list(c(1, 6.42)), 4)
  output <- capture.output(print(mb, overall = FALSE, details.methods = FALSE))
  expect_equal(sum(grepl("1, 6.42", output, fixed = TRUE)), 2)
})
