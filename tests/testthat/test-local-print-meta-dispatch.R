test_that("released legacy metabind objects are rejected before ordinary printing", {
  old <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), sm = "MD")
  old$version <- "8.2-1"
  class(old) <- c("metabind", "meta")
  for (printer in list(print, meta:::print.meta)) {
    error <- tryCatch(printer(old), error = identity)
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), "earlier than 9.0-0", fixed = TRUE)
    expect_match(conditionMessage(error), "Recreate it", fixed = TRUE)
    expect_match(conditionMessage(error), "Original call:", fixed = TRUE)
    expect_match(conditionMessage(error), "metagen(", fixed = TRUE)
  }
})

test_that("current metabind objects use the dedicated print method", {
  ma <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), sm = "MD")
  mb <- metabind(ma, update(ma, method.random.ci = "HK"))
  output <- capture.output(print(mb, details.methods = FALSE))
  dedicated <- capture.output(meta:::print.metabind(mb, details.methods = FALSE))
  expect_identical(output, dedicated)
  expect_true(any(output == "Meta-analyses results"))
})

test_that("ordinary single-study, excluded-study, and subgroup print paths still work", {
  models <- list(
    metabin(2, 20, 1, 20, method.tau = "DL"),
    metaprop(5, 50, method = "Inverse"),
    metabin(c(2, 4, 6), rep(20, 3), c(1, 3, 4), rep(20, 3),
              exclude = c(FALSE, TRUE, FALSE), method.tau = "DL"),
    metaprop(c(5, 10, 15, 20), rep(50, 4), method = "Inverse",
               subgroup = c("A", "A", "B", "B"), tau.common = TRUE))
  for (ma in models) {
    output <- capture.output(print(ma, details.methods = FALSE))
    direct <- capture.output(meta:::print.meta(ma, details.methods = FALSE))
    expect_identical(output, direct)
    expect_true(any(grepl("Number of observations:", output, fixed = TRUE)))
    expect_true(any(grepl("Number of events:", output, fixed = TRUE)))
  }
})
