test_that("grouped and individual print/digits options give identical output", {
  ma <- metagen(c(0.23456, 0.34567, 0.45678), c(0.12345, 0.23456, 0.34567),
                 sm = "MD", prediction = TRUE)
  mb <- metabind(ma, update(ma, method.random.ci = "HK", method.predict = "HK",
                            adhoc.hakn.pi = "se"))
  flags <- list(func = TRUE, method = TRUE, method.interval = TRUE,
                sm = TRUE, se = TRUE, pval = TRUE, stat = TRUE,
                level.ma = TRUE, level.predict = TRUE, level.hetstat = TRUE,
                tau2 = TRUE, tau2.ci = TRUE, tau = TRUE, tau.ci = TRUE,
                I2 = TRUE, I2.ci = TRUE, Q = TRUE, pval.Q = TRUE,
                H = TRUE, H.ci = TRUE)
  precision <- list(digits = 4, se = 1, stat = 3, pval = 1, pval.Q = 3,
                    Q = 1, tau2 = 3, tau = 2, I2 = 2)

  grouped <- capture.output(print(mb, print = flags, digits = precision,
                                  details.methods = FALSE))
  individual <- c(
    list(x = mb, digits = precision$digits, details.methods = FALSE),
    setNames(flags, paste0("print.", names(flags))),
    setNames(precision[-1], paste0("digits.", names(precision)[-1]))
  )
  legacy <- capture.output(do.call(print, individual))
  expect_identical(grouped, legacy)
  expect_true(any(grepl(formatC(ma$TE.common, digits = 4, format = "f"),
                        grouped, fixed = TRUE)))
})

test_that("automatic confidence-level selection stays independent by table", {
  ma <- metagen(1:3, 3:1, sm = "MD", common = TRUE, random = FALSE,
                 level.ma = 0.95)
  other <- metagen(1:3, 3:1, sm = "MD", common = FALSE, random = TRUE,
                    level.ma = 0.99)
  third <- metagen(1:3, 3:1, sm = "MD", common = FALSE, random = TRUE,
                    level.ma = 0.90)
  output <- capture.output(print(metabind(ma, other, third),
                                  details.methods = FALSE))
  common <- which(output == "Common effect model")
  random <- which(output == "Random effects model")
  heterogeneity <- which(output == "Information on heterogeneity")
  expect_false(any(grepl("level", output[common:(random - 1)], fixed = TRUE)))
  expect_true(any(grepl("level", output[random:(heterogeneity - 1)],
                       fixed = TRUE)))
})
