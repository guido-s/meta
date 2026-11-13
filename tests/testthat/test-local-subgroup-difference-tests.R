test_that("forest subgroup F labels retain both DF in all layouts", {
  ma <- subgroup_difference_fixture()
  for (layout in c("RevMan5", "meta", "BMJ", "JAMA")) {
    labels <- subgroup_difference_forest_label(ma, layout)
    expect_length(labels, 1)
    label <- labels[[1]]
    if (layout %in% c("RevMan5", "BMJ")) {
      expect_identical(label[[3]], "F")
      expect_match(label[[4]], "df = 1, 8", fixed = TRUE)
    }
    else {
      expect_identical(label[[3]][[2]], quote(F))
      expect_identical(label[[3]][[3]], "1, 8")
    }
    text <- paste(deparse(label), collapse = " ")
    expect_match(text, "4.86", fixed = TRUE)
    expect_match(text, "0.0587", fixed = TRUE)
  }
})

test_that("common-effect decomposition retains its within-group row", {
  ma <- update(subgroup_difference_fixture(c("classic", "HK"), method = "Inverse"),
                 common = TRUE)
  output <- capture.output(print(ma, details.methods = FALSE))
  tests <- grep("^Test for subgroup differences", output)
  expect_length(tests, 3)
  expect_equal(sum(grepl("^Within groups", output)), 1)
  expect_true(grepl("^Within groups", output[tests[1] + 3]))
  expect_false(any(grepl("^Within groups", output[seq.int(tests[2], length(output))])))
  expect_equal(sum(output == "Test of residual heterogeneity:"), 1)
})
