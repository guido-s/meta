test_that("prediction-only printing includes method details", {
  ma <- metagen(log(2:4), c(0.1, 0.2, 0.3), sm = "RR", prediction = TRUE)
  mb <- metabind(ma, update(ma, method.predict = "HK", adhoc.hakn.pi = "se"))

  output <- capture.output(print(mb, common = FALSE, random = FALSE))
  expect_true(any(grepl("Details of meta-analysis methods", output, fixed = TRUE)))
  expect_true(any(grepl("Hartung-Knapp", output, fixed = TRUE)))
  expect_false(any(grepl("- Inverse variance method", output, fixed = TRUE)))
})

test_that("prediction-only proportion printing retains component metadata", {
  ma <- metaprop(c(2, 4, 8), c(100, 200, 300), method = "Inverse",
                 prediction = TRUE, null.effect = 0.1)
  mb <- metabind(ma, update(ma, method.predict = "HK"))

  output <- capture.output(print(mb, common = FALSE, random = FALSE))
  expect_true(any(grepl("Meta-analysis of single proportions", output,
                       fixed = TRUE)))
})
