test_that("metacum selects random effects when prediction intervals are requested", {
  data(Fleiss1993bin)
  
  ma <- metabin(d.asp, n.asp, d.plac, n.plac,
                data = Fleiss1993bin, studlab = study, sm = "RR",
                method = "I", prediction = TRUE)
  
  mc <- metacum(ma)
  
  expect_identical(mc$pooled, "random")
  expect_false(mc$common)
  expect_true(mc$random)
  expect_true(mc$prediction)
})

test_that("metacum prediction argument selects random effects", {
  data(Fleiss1993bin)
  
  ma <- metabin(d.asp, n.asp, d.plac, n.plac,
                data = Fleiss1993bin, studlab = study, sm = "RR",
                method = "I", prediction = FALSE)
  
  mc <- metacum(ma, prediction = TRUE)
  
  expect_identical(mc$pooled, "random")
  expect_false(mc$common)
  expect_true(mc$random)
  expect_true(mc$prediction)
})

test_that("metacum explicit common effects selection is not overridden", {
  data(Fleiss1993bin)
  
  ma <- metabin(d.asp, n.asp, d.plac, n.plac,
                data = Fleiss1993bin, studlab = study, sm = "RR",
                method = "I", prediction = TRUE)
  
  mc <- metacum(ma, pooled = "common")
  
  expect_identical(mc$pooled, "common")
  expect_true(mc$common)
  expect_false(mc$random)
  expect_false(mc$prediction)
})

test_that("metacum keeps common effects default without prediction intervals", {
  data(Fleiss1993bin)
  
  ma <- metabin(d.asp, n.asp, d.plac, n.plac,
                data = Fleiss1993bin, studlab = study, sm = "RR",
                method = "I", prediction = FALSE)
  
  mc <- metacum(ma)
  
  expect_identical(mc$pooled, "common")
  expect_true(mc$common)
  expect_false(mc$random)
  expect_false(mc$prediction)
})
