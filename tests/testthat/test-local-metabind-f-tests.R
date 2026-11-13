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
