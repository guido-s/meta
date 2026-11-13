test_that("scalar tau confidence intervals use their own upper-bound signs", {
  x <- three_level_forest_fixture(three.level = FALSE)
  x$lower.tau2 <- 0.01
  x$upper.tau2 <- 0.81
  x$lower.tau <- 0.1
  x$upper.tau <- 0.9
  x$lower.tau2.resid <- 0.04
  x$upper.tau2.resid <- 0.64
  x$lower.tau.resid <- 0.2
  x$upper.tau.resid <- 0.8
  x$sign.lower.tau <- x$sign.lower.tau.resid <- ""
  x$sign.upper.tau <- ">"
  x$sign.upper.tau.resid <- "<"
  before <- x
  for (ci in c("tau2", "tau")) {
    fields <- paste0("hetstat.", ci, c("", ".resid"))
    drawn <- three_level_forest_capture(
      x, fields, resid.hetstat = TRUE, print = list(tau2 = TRUE, tau = TRUE,
        tau2.ci = ci == "tau2", tau.ci = ci == "tau"))
    expect_match(drawn$values[[fields[1]]],
                 if (ci == "tau2") "; >0.8100]" else "; >0.9000]", fixed = TRUE)
    expect_match(drawn$values[[fields[2]]],
                 if (ci == "tau2") "; <0.6400]" else "; <0.8000]", fixed = TRUE)
  }
  expect_identical(x, before)
})

test_that("merged forest placeholders handle hidden models and nested merges", {
  for (three.level in c(FALSE, TRUE)) {
    first <- three_level_forest_fixture(three.level = three.level)
    second <- three_level_forest_fixture(c("classic", "HK"), rho = 0.7,
                                          three.level = three.level)
    x <- metamerge(metamerge(first, second), first)
    before <- x
    for (selection in list(list(common = FALSE), list(random = FALSE),
                           list(common = FALSE, random = FALSE),
                           list(overall = FALSE), list(subgroup = c(TRUE, FALSE)))) {
      expect_silent(drawn <- do.call(
        three_level_forest_capture,
        c(list(x = x, fields = c("TEs", "lowTEs", "uppTEs")), selection)))
      expect_length(drawn$values$TEs, length(drawn$values$lowTEs))
      expect_length(drawn$values$TEs, length(drawn$values$uppTEs))
    }
    expect_identical(x, before)
  }
})
