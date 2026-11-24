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

test_that("metacum selects existing CI and prediction specifications", {
  ma <- metagen(c(-0.8, 0.9, -0.4, 1.2, 0.3, 1.6),
                c(0.12, 0.3, 0.18, 0.35, 0.22, 0.4),
                common = FALSE, prediction = TRUE, method.tau = "REML",
                method.random.ci = c("classic", "HK", "HK", "KR"),
                adhoc.hakn.ci = c("", "se"),
                method.predict = c("V", "S", "HK", "HK"),
                adhoc.hakn.pi = c("", "se"))
  original <- serialize(ma, NULL)
  fields <- c("TE", "seTE", "lower", "upper", "statistic", "pval",
              "df.random", "lower.predict", "upper.predict", "df.predict",
              "TE.pooled", "seTE.pooled", "lower.pooled", "upper.pooled",
              "statistic.pooled", "pval.pooled", "df.random.pooled",
              "lower.predict.pooled", "upper.predict.pooled",
              "df.predict.pooled", "text.pooled", "text.predict",
              "method.random.ci", "adhoc.hakn.ci",
              "method.predict", "adhoc.hakn.pi")

  for (i in seq_along(ma$method.random.ci)) {
    for (j in seq_along(ma$method.predict)) {
      method.ci <- unname(ma$method.random.ci[i])
      adhoc.ci <- unname(ma$adhoc.hakn.ci[i])
      method.pi <- unname(ma$method.predict[j])
      adhoc.pi <- unname(ma$adhoc.hakn.pi[j])

      selected <- metacum(ma,
                          method.random.ci = method.ci,
                          adhoc.hakn.ci = adhoc.ci,
                          method.predict = method.pi,
                          adhoc.hakn.pi = adhoc.pi)
      reference <- metacum(update(ma,
                                   method.random.ci = method.ci,
                                   adhoc.hakn.ci = adhoc.ci,
                                   method.predict = method.pi,
                                   adhoc.hakn.pi = adhoc.pi,
                                   text.random = unname(ma$text.random[i]),
                                   text.predict = unname(ma$text.predict[j])))

      for (field in fields)
        expect_equal(unname(selected[[field]]), unname(reference[[field]]),
                     info = paste(method.ci, adhoc.ci, method.pi, adhoc.pi,
                                  field))
    }
  }
  expect_identical(serialize(ma, NULL), original)
})

test_that("metacum retains the shared GLMM estimate when selecting a CI method", {
  ma <- metaprop(c(10, 20, 30, 40, 50), rep(100, 5),
                 method = "GLMM", common = FALSE,
                 method.random.ci = c("classic", "HK"))
  expect_length(ma$TE.random, 1L)

  for (i in seq_along(ma$method.random.ci)) {
    selected <- metacum(ma, method.random.ci = ma$method.random.ci[i])
    expect_equal(unname(selected$TE.pooled), unname(ma$TE.random))
    expect_equal(unname(selected$seTE.pooled), unname(ma$seTE.random[i]))
    expect_equal(unname(selected$lower.pooled), unname(ma$lower.random[i]))
    expect_equal(unname(selected$upper.pooled), unname(ma$upper.random[i]))
  }
})

test_that("metacum defaults to the first available interval settings", {
  ma <- metagen(c(-0.8, 0.9, -0.4, 1.2, 0.3, 1.6),
                c(0.12, 0.3, 0.18, 0.35, 0.22, 0.4),
                common = FALSE, prediction = TRUE,
                method.random.ci = c("classic", "HK", "HK"),
                adhoc.hakn.ci = c("se", ""),
                method.predict = c("V", "HK", "HK"),
                adhoc.hakn.pi = c("se", ""))

  default <- metacum(ma)
  explicit <- metacum(ma, method.random.ci = "classic", method.predict = "V")
  fields <- setdiff(names(default), "call")
  expect_equal(default[fields], explicit[fields])

  selected <- metacum(ma, method.random.ci = "HK", method.predict = "HK")
  expect_equal(unname(selected$adhoc.hakn.ci), "se")
  expect_equal(unname(selected$adhoc.hakn.pi), "se")

  common <- metacum(ma, pooled = "common")
  expect_true(is.na(common$df.random.pooled))
  expect_true(is.na(common$df.predict.pooled))
})

test_that("metacum validates interval selectors and warns about ignored no", {
  ma <- metagen(c(-0.8, 0.9, -0.4, 1.2, 0.3, 1.6),
                c(0.12, 0.3, 0.18, 0.35, 0.22, 0.4),
                common = FALSE, prediction = TRUE,
                method.random.ci = c("classic", "HK"),
                method.predict = c("V", "HK"))

  expect_error(metacum(ma, method.random.ci = "KR"), "method.random.ci")
  expect_error(metacum(ma, method.predict = "S"), "method.predict")
  expect_error(metacum(ma, method.random.ci = "classic", adhoc.hakn.ci = "se"),
               "adhoc.hakn.ci")
  expect_error(metacum(ma, method.predict = "V", adhoc.hakn.pi = "se"),
               "adhoc.hakn.pi")

  expect_warning(ignored <- metacum(ma, no = 2), "Argument 'no' is ignored")
  expect_equal(ignored$lower, metacum(ma)$lower)
  expect_equal(ignored$lower.predict, metacum(ma)$lower.predict)
})
