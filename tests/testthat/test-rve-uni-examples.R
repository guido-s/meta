test_that("RVE methods work in univariate random effects meta-analysis", {
  data(Fleiss1993bin)
  data(Fleiss1993cont)
  data(smoking)

  mbin <- metabin(d.asp, n.asp, d.plac, n.plac, study,
                  data = Fleiss1993bin,
                  method = "Inverse", method.random.ci = "CR1",
                  random = TRUE, common = FALSE)
  mcont <- metacont(n.psyc, mean.psyc, sd.psyc,
                    n.cont, mean.cont, sd.cont, study,
                    data = Fleiss1993cont,
                    method.random.ci = "CR1",
                    random = TRUE, common = FALSE)
  mcor <- metacor(c(0.1, 0.2, 0.3, 0.4), c(20, 30, 40, 50),
                  studlab = paste("study", 1:4),
                  method.random.ci = "CR1",
                  random = TRUE, common = FALSE)
  mgen <- metagen(c(1, 2, 3, 4), c(0.2, 0.3, 0.4, 0.5),
                  studlab = paste("study", 1:4),
                  method.random.ci = "CR1",
                  random = TRUE, common = FALSE)
  minc <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
                  study, data = smoking,
                  method = "Inverse", method.random.ci = "CR1",
                  random = TRUE, common = FALSE)
  mmean <- metamean(c(20, 30, 40, 50), c(1, 2, 3, 4), c(0.2, 0.3, 0.4, 0.5),
                    studlab = paste("study", 1:4),
                    method.random.ci = "CR1",
                    random = TRUE, common = FALSE)
  mprop <- metaprop(c(1, 5, 10, 15), c(20, 30, 40, 50),
                    studlab = paste("study", 1:4),
                    method = "Inverse", method.random.ci = "CR1",
                    random = TRUE, common = FALSE)
  mrate <- metarate(c(1, 5, 10, 15), c(20, 30, 40, 50),
                    studlab = paste("study", 1:4),
                    method = "Inverse", method.random.ci = "CR1",
                    random = TRUE, common = FALSE)

  res <- list(mbin, mcont, mcor, mgen, minc, mmean, mprop, mrate)

  expect_true(all(vapply(res, function(x) x$method.random.ci == "CR1", TRUE)))
  expect_true(all(vapply(res, function(x) is.finite(x$seTE.random), TRUE)))
  expect_true(all(vapply(res, function(x) is.finite(x$df.random), TRUE)))
})

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
