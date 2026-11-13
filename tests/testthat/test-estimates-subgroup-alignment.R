estimates_subgroup_fixture <- function(methods, backend = "Inverse") {
  if (backend == "GLMM") {
    dat <- data.frame(event = c(10, 20, 30, 40, 50, 30, 40, 50, 60, 70),
                      total = 100, group = rep(c("A", "B"), each = 5))
    metaprop(event, total, data = dat, method = "GLMM", common = TRUE,
               subgroup = group, tau.common = TRUE, method.random.ci = methods)
  }
  else
    metagen(c(-0.6, -0.2, 0.15, 0.4, 0.8, 1.1, 1.4, 1.9),
              c(0.12, 0.18, 0.15, 0.2, 0.14, 0.22, 0.19, 0.16),
              subgroup = rep(c("A", "B"), c(3, 5)), sm = "MD",
              method.tau = "DL", method.random.ci = methods)
}

expect_estimates_subgroup_alignment <- function(x, model, extracted = NULL) {
  if (is.null(extracted))
    extracted <- estimates(
      x, study.results = FALSE, overall = FALSE,
      common = model == "common", random = model == "random",
      prediction = FALSE, prediction.subgroup = FALSE,
      se = TRUE, statistic = TRUE, pval = TRUE, backtransf = FALSE,
      digits = 12, digits.se = 12, digits.stat = 12, digits.pval = 12,
      writexl = FALSE)
  labels <- x[[paste0("text.", model)]]
  methods <- x[[paste0("method.", model, ".ci")]]
  expect_length(labels, length(methods))
  expect_equal(nrow(extracted), length(x$subgroup.levels) * length(labels))
  expect_equal(sum(names(extracted) == "estimate"), 1)
  expect_false(any(grepl("^estimate[.]", names(extracted))))
  expect_identical(attr(extracted, "x")[[paste0("method.", model, ".ci")]], methods)
  fields <- c(estimate = "TE", se = "seTE", lower = "lower", upper = "upper",
               statistic = "statistic", pval = "pval")

  for (g in seq_along(x$subgroup.levels)) {
    for (j in seq_along(labels)) {
      row <- (g - 1L) * length(labels) + j
      info <- paste(class(x)[1], model, x$subgroup.levels[g], methods[j])
      expect_identical(unname(extracted$subgroup[row]), unname(x$subgroup.levels[g]),
                       info = info)
      expect_identical(unname(extracted$studlab[row]), unname(labels[j]), info = info)
      expect_equal(unname(extracted$k[row]), unname(x$k.w[g]), info = info)
      for (field in names(fields)) {
        stored <- x[[paste0(fields[[field]], ".", model, ".w")]]
        expected <- if (is.matrix(stored)) stored[g, j] else stored[g]
        expect_equal(unname(extracted[[field]][row]),
                      unname(round(expected, 12)), tolerance = 1e-12, info = info)
      }
      if (model == "random" && "df" %in% names(extracted)) {
        stored <- x$df.random.w
        expected <- if (is.matrix(stored)) stored[g, j] else stored[g]
        expect_equal(unname(extracted$df[row]), unname(round(expected, 2)), info = info)
      }
    }
  }
  invisible(extracted)
}

test_that("random subgroup extraction aligns all fields for single and reordered methods", {
  for (backend in c("Inverse", "GLMM")) {
    for (methods in list("classic", "HK", c("classic", "HK"), c("HK", "classic"))) {
      ma <- estimates_subgroup_fixture(methods, backend)
      expect_gt(diff(range(ma$TE.random.w)), 0.1)
      expect_estimates_subgroup_alignment(ma, "random")
    }
  }
})

test_that("common subgroup extraction aligns both supported single CI methods", {
  for (method in c("classic", "IVhet")) {
    ma <- update(estimates_subgroup_fixture("classic"), method.common.ci = method)
    expect_estimates_subgroup_alignment(ma, "common")
  }
  expect_estimates_subgroup_alignment(estimates_subgroup_fixture("HK", "GLMM"),
                                      "common")
})

test_that("subgroup extraction handles vector and equivalent matrix estimate storage", {
  ma <- estimates_subgroup_fixture(c("classic", "HK"))
  for (model in c("common", "random")) {
    vector <- expect_estimates_subgroup_alignment(ma, model)
    matrix.object <- ma
    field <- paste0("TE.", model, ".w")
    labels <- ma[[paste0("text.", model)]]
    matrix.object[[field]] <- matrix(
      rep(ma[[field]], length(labels)), nrow = length(ma$subgroup.levels),
      dimnames = list(ma$subgroup.levels, ma[[paste0("method.", model, ".ci")]]))
    expanded <- expect_estimates_subgroup_alignment(matrix.object, model)
    for (column in names(vector))
      expect_identical(expanded[[column]], vector[[column]])
  }
})

test_that("combined common and random extraction preserves model blocks and method order", {
  for (methods in list(c("classic", "HK"), c("HK", "classic"))) {
    ma <- estimates_subgroup_fixture(methods)
    extracted <- estimates(
      ma, study.results = FALSE, overall = FALSE, common = TRUE, random = TRUE,
      prediction = FALSE, prediction.subgroup = FALSE,
      se = TRUE, statistic = TRUE, pval = TRUE, backtransf = FALSE,
      digits = 12, digits.se = 12, digits.stat = 12, digits.pval = 12,
      writexl = FALSE)
    n.common <- length(ma$subgroup.levels) * length(ma$text.common)
    expect_estimates_subgroup_alignment(ma, "common", extracted[seq_len(n.common), ])
    expect_estimates_subgroup_alignment(ma, "random",
                                        extracted[-seq_len(n.common), ])
  }
})
