test_that("IRFT restrictions avoid ambiguous Freeman-Tukey times", {
  m <- metarate(c(2, 5, 8, 12), c(30, 40, 35, 50), sm = "IRFT")
  
  expect_error(metaadd(m, type = "random", TE = 0.2, lower = 0.1, upper = 0.3),
               "not available")
  expect_error(metagen(c(0.1, 0.2), c(0.01, 0.02), sm = "IRFT",
                       transf = FALSE),
               "not available")
})
