test_that("metainf selects random effects when prediction intervals are requested", {
  data(Fleiss1993bin)
  
  ma <- metabin(d.asp, n.asp, d.plac, n.plac,
                data = Fleiss1993bin, studlab = study, sm = "RR",
                method = "I", prediction = TRUE)
  
  mi <- metainf(ma)
  
  expect_identical(mi$pooled, "random")
  expect_false(mi$common)
  expect_true(mi$random)
  expect_true(mi$prediction)
})

test_that("metainf prediction argument selects random effects", {
  data(Fleiss1993bin)
  
  ma <- metabin(d.asp, n.asp, d.plac, n.plac,
                data = Fleiss1993bin, studlab = study, sm = "RR",
                method = "I", prediction = FALSE)
  
  mi <- metainf(ma, prediction = TRUE)
  
  expect_identical(mi$pooled, "random")
  expect_false(mi$common)
  expect_true(mi$random)
  expect_true(mi$prediction)
})

test_that("metainf explicit common effects selection is not overridden", {
  data(Fleiss1993bin)
  
  ma <- metabin(d.asp, n.asp, d.plac, n.plac,
                data = Fleiss1993bin, studlab = study, sm = "RR",
                method = "I", prediction = TRUE)
  
  mi <- metainf(ma, pooled = "common")
  
  expect_identical(mi$pooled, "common")
  expect_true(mi$common)
  expect_false(mi$random)
  expect_false(mi$prediction)
})

test_that("metainf keeps common effects default without prediction intervals", {
  data(Fleiss1993bin)
  
  ma <- metabin(d.asp, n.asp, d.plac, n.plac,
                data = Fleiss1993bin, studlab = study, sm = "RR",
                method = "I", prediction = FALSE)
  
  mi <- metainf(ma)
  
  expect_identical(mi$pooled, "common")
  expect_true(mi$common)
  expect_false(mi$random)
  expect_false(mi$prediction)
})
