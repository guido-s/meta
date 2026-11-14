test_that("metainf prediction selection does not depend on input random flag", {
  for (inherited in c(FALSE, TRUE)) {
    ma <- metagen(c(-0.8, 0.9, -0.4, 1.2, 0.3, 1.6),
                  c(0.12, 0.3, 0.18, 0.35, 0.22, 0.4),
                  common = TRUE, random = FALSE, prediction = inherited,
                  method.tau = "DL", method.random.ci = "classic",
                  method.predict = "V")
    original <- serialize(ma, NULL)
    expect_true(ma$common)
    expect_false(ma$random)
    expect_identical(ma$prediction, inherited)

    requested <- metainf(ma, prediction = TRUE)
    default <- metainf(ma)
    disabled <- metainf(ma, prediction = FALSE)
    common <- metainf(ma, pooled = "common")
    random <- metainf(ma, pooled = "random")
    explicit.common <- metainf(ma, pooled = "common", prediction = TRUE)
    explicit.random <- metainf(ma, pooled = "random", prediction = FALSE)

    expect_identical(requested$pooled, "random")
    expect_false(requested$common)
    expect_true(requested$random)
    expect_true(requested$prediction)
    expect_identical(default$pooled, if (inherited) "random" else "common")
    expect_identical(default$prediction, inherited)
    expect_identical(disabled$pooled, "common")
    expect_true(disabled$common)
    expect_false(disabled$random)
    expect_false(disabled$prediction)
    expect_identical(common$pooled, "common")
    expect_true(common$common)
    expect_false(common$random)
    expect_false(common$prediction)
    expect_identical(random$pooled, "random")
    expect_false(random$common)
    expect_true(random$random)
    expect_identical(random$prediction, inherited)
    expect_identical(explicit.common$pooled, "common")
    expect_identical(explicit.random$pooled, "random")
    expect_false(explicit.random$prediction)
    expect_identical(serialize(ma, NULL), original)
  }
})

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
