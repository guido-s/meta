test_that("logical-list validation supports subgroup vectors and optional NULLs", {
  flags <- list(main = TRUE, subgroup = c(TRUE, FALSE), optional = NULL)
  expect_silent(meta:::chklogical_list(flags, NULL.ok = "optional",
                                      length = list(subgroup = 0)))
  expect_error(meta:::chklogical_list(flags, NULL.ok = "optional"), "flags[.]subgroup")
  flags$subgroup[2] <- NA
  expect_error(meta:::chklogical_list(flags, NULL.ok = "optional",
                                      length = list(subgroup = 0)), "flags[.]subgroup")
})

test_that("print lists retain legacy display settings and vector CI controls", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3), method.tau = "DL")
  settings <- list(I2 = TRUE, I2.ci = TRUE, tau2 = TRUE, tau = FALSE,
                    Q = TRUE, pval.Q = TRUE, Rb = TRUE, stat = TRUE,
                    subgroup.name = FALSE, tau2.ci.subgroup = c(TRUE, FALSE),
                    tau.ci.subgroup = c(FALSE, TRUE))
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      observed <<- list(...)
      do.call(draw, observed)
    }, .package = "meta")
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())

  baseline <- NULL
  for (legacy in c(FALSE, TRUE)) {
    if (legacy)
      options <- setNames(settings, paste0("print.", names(settings)))
    else
      options <- list(print = settings)
    do.call(meta::forest, c(list(x = ma), options))
    labels <- lapply(observed$cols$col.studlab$labels, function(grob)
      paste(deparse(grob$label), collapse = " "))
    if (!legacy)
      baseline <- labels
    else
      expect_identical(labels, baseline)
  }
  expect_silent(forest(ma, print = list(I2 = NULL, tau2 = NULL, pval.Q = NULL,
                                        subgroup.name = NULL,
                                        tau2.ci.subgroup = NULL, tau.ci.subgroup = NULL)))
  expect_warning(forest(ma, print.zval = FALSE), "print.stat")
  expect_error(forest(ma, print = list(I2 = "invalid")), "print[.]I2")
  expect_error(forest(ma, print = list(I2 = NA)), "print[.]I2")
  expect_error(forest(ma, print = list(stat = NULL)), "print[.]stat")
  expect_error(forest(ma, print = list(Q = NULL)), "print[.]Q")
  expect_error(forest(ma, print = list(tau2.ci.subgroup = c(TRUE, NA))),
               "print[.]tau2.ci.subgroup")
  expect_error(forest(ma, print = list(tau2.ci.subgroup = rep(TRUE, 3))),
               "number of subgroups")
})
