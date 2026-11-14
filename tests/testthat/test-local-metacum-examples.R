test_that("metacum prediction selection does not depend on input random flag", {
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

    requested <- metacum(ma, prediction = TRUE)
    default <- metacum(ma)
    disabled <- metacum(ma, prediction = FALSE)
    common <- metacum(ma, pooled = "common")
    random <- metacum(ma, pooled = "random")
    explicit.common <- metacum(ma, pooled = "common", prediction = TRUE)
    explicit.random <- metacum(ma, pooled = "random", prediction = FALSE)

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
