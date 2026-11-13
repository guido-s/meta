subgroup_effect_fixture <- function(methods = c("classic", "HK")) {
  x <- metagen(seq(-0.3, 0.8, length.out = 12), rep(0.2, 12),
               subgroup = rep(c("A", "B", "C"), c(3, 4, 5)),
               method.tau = "DL", method.random.ci = methods)
  x$pval.common.w[] <- c(0.1111, 0.2222, 0.3333)
  x$statistic.common.w[] <- c(1.111, 2.222, 3.333)
  x$pval.random.w[] <- seq(0.1234, by = 0.1,
                           length.out = length(x$pval.random.w))
  x$statistic.random.w[] <- seq(11.111, by = 1.111,
                                length.out = length(x$statistic.random.w))
  x
}

subgroup_effect_labels <- function(x, ...) {
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      args <- list(...)
      observed <<- lapply(args$cols$col.studlab$labels, `[[`, "label")
      do.call(draw, args)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 16, height = 16)
  on.exit(grDevices::dev.off(), add = TRUE)
  forest(x, layout = "meta", test.effect.subgroup = TRUE,
         label.test.effect.subgroup.common = "COMMON_EFFECT",
         label.test.effect.subgroup.random = "RANDOM_EFFECT",
         digits.stat = 3, digits.pval = 4, scientific.pval = FALSE,
         zero.pval = TRUE, ...)
  lapply(c(common = "COMMON_EFFECT", random = "RANDOM_EFFECT"), function(tag)
    Filter(function(z) is.call(z) && identical(z[[2]], tag), observed))
}

expect_subgroup_effect <- function(label, pval, statistic, df = NULL) {
  expect_equal(unname(label[[5]]),
               unname(meta:::formatN(statistic, 3, "NA", big.mark = "")))
  expect_equal(trimws(unname(label[[8]])),
               trimws(unname(meta:::formatPT(pval, lab = TRUE, labval = "",
                                     digits = 4, zero = TRUE, JAMA = FALSE,
                                     scientific = FALSE, lab.NA = "NA"))))
  if (is.null(df))
    expect_equal(label[[3]], quote(italic(z)))
  else {
    expect_equal(label[[3]][[1]], as.name("["))
    expect_equal(label[[3]][[2]], quote(italic(t)))
    expect_equal(unname(label[[3]][[3]]), unname(df))
  }
}

test_that("multiple random methods use the first slot within each group", {
  for (methods in list(c("classic", "HK"), c("HK", "classic"))) {
    x <- subgroup_effect_fixture(methods)
    labels <- subgroup_effect_labels(x)
    expect_length(labels$common, 3)
    expect_length(labels$random, 3)
    for (i in 1:3) {
      expect_subgroup_effect(labels$common[[i]], x$pval.common.w[i],
                             x$statistic.common.w[i])
      expect_subgroup_effect(labels$random[[i]], x$pval.random.w[i, 1],
                             x$statistic.random.w[i, 1],
                             if (methods[1] == "HK") x$df.random.w[i, 1])
    }
  }
})

test_that("sparse enabled slots select their own method and degrees of freedom", {
  x <- subgroup_effect_fixture()
  x$df.random.w[, 2] <- c(2.25, 3.5, 4.75)
  enabled <- matrix(c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE), 3, byrow = TRUE)
  x$pval.random.w[2, 1] <- NA_real_
  # A missing first p-value and a disabled first slot both select HK.
  labels <- subgroup_effect_labels(x, random.subgroup = enabled)
  expect_length(labels$random, 3)
  for (i in 1:3)
    expect_subgroup_effect(labels$random[[i]], x$pval.random.w[i, 2],
                           x$statistic.random.w[i, 2], x$df.random.w[i, 2])

  # Both methods remain present, but different components enable different slots.
  enabled[2, ] <- c(FALSE, TRUE)
  enabled[3, ] <- c(TRUE, FALSE)
  labels <- subgroup_effect_labels(x, random.subgroup = enabled)
  for (i in 1:3) {
    slot <- c(2, 2, 1)[i]
    expect_subgroup_effect(labels$random[[i]], x$pval.random.w[i, slot],
                           x$statistic.random.w[i, slot],
                           if (slot == 2) x$df.random.w[i, slot])
  }
})

test_that("single-method subgroup test text is unchanged", {
  for (method in c("classic", "HK")) {
    x <- subgroup_effect_fixture(method)
    labels <- subgroup_effect_labels(x)
    expect_length(labels$random, 3)
    for (i in 1:3)
      expect_subgroup_effect(labels$random[[i]], x$pval.random.w[i],
                             x$statistic.random.w[i],
                             if (method == "HK") x$df.random.w[i])
  }
})

test_that("p-value-only rows use the selected random slot", {
  x <- subgroup_effect_fixture()
  enabled <- matrix(c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE), 3, byrow = TRUE)
  labels <- subgroup_effect_labels(x, random.subgroup = enabled,
                                   print.stat = FALSE)
  expect_length(labels$random, 3)
  for (i in 1:3) {
    pval <- x$pval.random.w[i, c(2, 1, 2)[i]]
    expect_equal(trimws(unname(labels$random[[i]][[5]])),
                 trimws(meta:::formatPT(pval, lab = TRUE, labval = "",
                                       digits = 4, zero = TRUE, JAMA = FALSE,
                                       scientific = FALSE, lab.NA = "NA")))
  }
})
