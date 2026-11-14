test_that("named transformation arguments remain intact and do not alter calculations", {
  fun <- ".meta_transformation_columns_affine"
  assign(fun, function(x, offset = 0, scale = 1) (x + offset) * scale,
         envir = .GlobalEnv)
  on.exit(rm(list = fun, envir = .GlobalEnv))
  args <- list(offset = 1, scale = 100)
  ma <- metagen_transformation_columns()
  custom <- update(ma, func.backtransf = fun, args.backtransf = args)
  mb <- metabind(custom, ma)
  original <- mb$overall
  expect_identical(mb$overall$args.backtransf, list(args, NULL))
  expect_identical(mb$overall$args.transf, list(NULL, NULL))
  expect_identical(meta:::format_args(c(args, list(unused = NULL))),
                   "offset = 1, scale = 100, unused = NULL")
  output <- capture.output(print(mb, details.methods = FALSE))
  expect_equal(sum(grepl("backtransf +args", output)), 3L)
  expect_equal(sum(grepl("offset = 1, scale = 100", output, fixed = TRUE)), 3L)
  shown <- forest_transformation_columns(mb, layout = "meta")
  hidden <- forest_transformation_columns(
    mb, layout = "meta", print = list(backtransf = FALSE, args.backtransf = TRUE))
  numeric <- c("TE", "lower", "upper")
  for (col in numeric) {
    expected <- shown$data[[paste0(col, ".orig")]]
    expected[shown$data$id == 1] <- (expected[shown$data$id == 1] + 1) * 100
    expect_equal(shown$data[[col]], expected)
  }
  expect_identical(shown$data[numeric], hidden$data[numeric])
  expect_equal(shown$data$args.backtransf,
               ifelse(shown$data$id == 1, "offset = 1, scale = 100", "."))
  expect_identical(mb$overall, original)
})

test_that("MD exponential and default rows retain their numeric plot data", {
  ma <- metagen_transformation_columns()
  mb <- metabind(update(ma, func.backtransf = "exp"), ma)
  res <- forest_transformation_columns(mb, layout = "meta")
  baseline <- forest_transformation_columns(metabind(ma, ma), layout = "meta")
  for (col in c("TE", "lower", "upper")) {
    expected <- baseline$data[[col]]
    expected[res$data$id == 1] <- exp(expected[res$data$id == 1])
    expect_equal(res$data[[col]], expected)
  }
  output <- capture.output(print(mb, digits = 2, details.methods = FALSE))
  expect_true(any(grepl(formatC(exp(ma$TE.common), digits = 2, format = "f"),
                        output, fixed = TRUE)))
  expect_identical(mb$overall$args.backtransf, list(NULL, NULL))
})
