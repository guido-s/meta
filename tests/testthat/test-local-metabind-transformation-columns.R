test_that("print automatically shows functions but not absent arguments", {
  ma <- metagen_transformation_columns()
  mb <- metabind(update(ma, func.backtransf = "exp"), ma)
  output <- capture.output(print(mb, print = list(k = FALSE), details.methods = FALSE))
  automatic <- capture.output(print(
    mb, print = list(backtransf = NULL, args.backtransf = NULL, k = FALSE),
    details.methods = FALSE))
  expect_identical(output, automatic)
  expect_equal(sum(grepl("backtransf", output, fixed = TRUE)), 3L)
  expect_false(any(grepl("args", output, fixed = TRUE)))
  expect_equal(sum(grepl("^meta2 +\\.", output)), 3L)
  het <- which(output == "Information on heterogeneity")
  expect_false(any(grepl("backtransf|args", output[het:length(output)])))
  plain <- capture.output(print(metabind(ma, ma), details.methods = FALSE))
  expect_false(any(grepl("backtransf|args", plain)))
})

test_that("print flags force or suppress columns and respect selected rows", {
  ma <- metagen_transformation_columns()
  mb <- metabind(update(ma, func.backtransf = "exp"), ma)
  flags <- list(list(backtransf = FALSE),
                list(args.backtransf = TRUE),
                list(backtransf = FALSE, args.backtransf = TRUE),
                list(backtransf = TRUE, args.backtransf = FALSE))
  headings <- vapply(flags, function(flag) {
    output <- capture.output(print(mb, print = flag, details.methods = FALSE))
    c(func = sum(grepl("backtransf", output, fixed = TRUE)),
      args = sum(grepl("args", output, fixed = TRUE)))
  }, c(func = 0L, args = 0L))
  expect_equal(unname(headings), rbind(c(0L, 3L, 0L, 3L), c(0L, 3L, 0L, 0L)))
  forced <- capture.output(print(
    metabind(ma, ma), print = list(backtransf = TRUE, args.backtransf = TRUE, k = FALSE),
    details.methods = FALSE))
  expect_equal(sum(grepl("backtransf +args", forced)), 3L)
  expect_equal(sum(grepl("^meta[12] +\\. +\\.", forced)), 6L)

  selected <- metabind(update(ma, common = FALSE, func.backtransf = "exp"), ma)
  output <- capture.output(print(selected, details.methods = FALSE))
  expect_equal(sum(grepl("backtransf", output, fixed = TRUE)), 2L)
  hidden <- capture.output(print(selected, random = FALSE, prediction = FALSE,
                                 details.methods = FALSE))
  expect_false(any(grepl("backtransf|args", hidden)))
})

test_that("forest automatic columns and explicit flags use selected rows", {
  ma <- metagen_transformation_columns()
  mb <- metabind(update(ma, func.backtransf = "exp"), ma)
  cols <- c("func.backtransf", "args.backtransf")
  flags <- list(NULL, list(backtransf = NULL, args.backtransf = NULL),
                list(backtransf = FALSE), list(args.backtransf = TRUE),
                list(backtransf = FALSE, args.backtransf = TRUE))
  present <- vapply(flags, function(flag) {
    res <- forest_transformation_columns(mb, layout = "meta", print = flag)
    cols %in% res$leftcols
  }, logical(2))
  expect_equal(present, rbind(c(TRUE, TRUE, FALSE, TRUE, FALSE),
                              c(FALSE, FALSE, FALSE, TRUE, FALSE)))
  plain <- metabind(ma, ma)
  expect_false(any(cols %in%
                     forest_transformation_columns(plain, layout = "meta")$leftcols))
  forced <- forest_transformation_columns(
    plain, layout = "meta", print = list(backtransf = TRUE, args.backtransf = TRUE))
  expect_true(all(cols %in% forced$leftcols))
  expect_true(all(forced$data[cols] == "."))
  selected <- metabind(update(ma, common = FALSE, func.backtransf = "exp"), ma)
  res <- forest_transformation_columns(selected, layout = "meta",
                                        random = FALSE, prediction = FALSE)
  expect_false(any(cols %in% res$leftcols))
})

test_that("forest layouts order transformation columns and preserve explicit labels", {
  fun <- ".meta_transformation_columns_layout_affine"
  assign(fun, function(x, offset = 0, scale = 1) (x + offset) * scale,
         envir = .GlobalEnv)
  on.exit(rm(list = fun, envir = .GlobalEnv))
  ma <- metagen_transformation_columns()
  mb <- metabind(update(ma, func.backtransf = fun,
                        args.backtransf = list(offset = 1, scale = 100)), ma)
  cols <- c("func.backtransf", "args.backtransf")
  for (layout in c("meta", "meta8", "BMJ", "RevMan5", "JAMA")) {
    res <- forest_transformation_columns(mb, layout = layout)
    effects <- if (layout %in% c("BMJ", "RevMan5", "JAMA"))
      "effect.ci" else c("effect", "ci")
    right <- layout %in% c("meta8", "BMJ")
    displayed <- if (right) res$rightcols else res$leftcols
    expected <- if (right) c(effects, cols) else c(cols, effects)
    expect_identical(displayed[displayed %in% c(cols, effects)], expected,
                     info = layout)
  }
  expect_equal(res$leftlabs[match(cols, res$leftcols)],
               c("Back-\ntransformation", "Arguments"))
  res <- forest_transformation_columns(
    mb, layout = "meta", leftcols = c("studlab", "args.backtransf"),
    leftlabs = c("Analysis", "Parameters"),
    rightcols = c("func.backtransf", "effect"),
    rightlabs = c("Transform", "Estimate"), print = list(backtransf = FALSE))
  expect_identical(list(res$leftcols, res$leftlabs, res$rightcols, res$rightlabs),
                   list(c("studlab", "args.backtransf"), c("Analysis", "Parameters"),
                        c("func.backtransf", "effect"), c("Transform", "Estimate")))
})
