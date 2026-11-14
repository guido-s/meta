test_that("RVE methods are rejected for unsupported random effects methods", {
  data(Fleiss1993bin)
  data(smoking)

  expect_error(
    metabin(d.asp, n.asp, d.plac, n.plac, study,
            data = Fleiss1993bin, sm = "OR",
            method = "Peto", method.random.ci = "CR1"),
    "not available"
  )
  expect_error(
    metabin(d.asp, n.asp, d.plac, n.plac, study,
            data = Fleiss1993bin, sm = "OR",
            method = "GLMM", method.random.ci = "CR1"),
    "not available"
  )
  expect_error(
    metabin(d.asp, n.asp, d.plac, n.plac, study,
            data = Fleiss1993bin, sm = "OR",
            method = "LRP", method.random.ci = "CR1"),
    "not available"
  )
  expect_error(
    metabin(d.asp, n.asp, d.plac, n.plac, study,
            data = Fleiss1993bin, sm = "OR",
            method = "SSW", method.random.ci = "CR1"),
    "not available"
  )
  expect_error(
    metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
            study, data = smoking,
            method = "GLMM", method.random.ci = "CR1"),
    "not available"
  )
  expect_error(
    metaprop(c(1, 5, 10, 15), c(20, 30, 40, 50),
             method = "GLMM", method.random.ci = "CR1"),
    "not available"
  )
  expect_error(
    metarate(c(1, 5, 10, 15), c(20, 30, 40, 50),
             method = "GLMM", method.random.ci = "CR1"),
    "not available"
  )
})
