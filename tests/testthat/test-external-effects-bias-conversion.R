test_that("OR and MD study effects agree with independent escalc oracles", {
  # escalc is imported but not actively called by meta: these are numerical
  # oracles, not coverage of an internal escalc call.
  event.e <- c(8, 15, 21)
  event.c <- c(12, 11, 18)
  n.e <- c(45, 62, 78)
  n.c <- c(51, 58, 83)
  binary <- metabin(event.e, n.e, event.c, n.c,
                    sm = "OR", method = "Inverse", method.tau = "DL")
  ref.or <- metafor::escalc("OR", ai = event.e, bi = n.e - event.e,
                            ci = event.c, di = n.c - event.c)
  expect_equal(binary$TE, as.numeric(ref.or$yi))
  expect_equal(binary$seTE^2, as.numeric(ref.or$vi))

  mean.e <- c(4.2, 6.5, 8.1)
  mean.c <- c(3.7, 7.2, 6.8)
  sd.e <- c(1.1, 2.3, 1.7)
  sd.c <- c(1.8, 1.5, 2.6)
  continuous <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
                         sm = "MD", pooledvar = FALSE, method.tau = "DL")
  ref.md <- metafor::escalc("MD", n1i = n.e, m1i = mean.e, sd1i = sd.e,
                            n2i = n.c, m2i = mean.c, sd2i = sd.c)
  expect_equal(continuous$TE, as.numeric(ref.md$yi))
  expect_equal(continuous$seTE^2, as.numeric(ref.md$vi))
})

test_that("metabias Egger follows the metafor lm regression reference", {
  TE <- c(-.7, -.2, .1, .6, 1, 1.5)
  seTE <- c(.15, .22, .18, .30, .20, .25)
  ma <- metagen(TE, seTE, method.tau = "DL")
  bias <- metabias(ma, method.bias = "Egger", k.min = 3)
  base <- metafor::rma.uni(yi = TE, sei = seTE, method = "DL", test = "t")
  ref <- metafor::regtest(base, predictor = "sei", model = "lm")
  coefficients <- coef(summary(ref$fit))

  expect_equal(unname(bias$statistic), as.numeric(ref$zval))
  expect_equal(bias$df, ref$dfs)
  expect_equal(bias$pval, ref$pval)
  expect_equal(bias$p.value, ref$pval)
  expect_equal(unname(bias$estimate), unname(coefficients[2, 1:2]))
  expect_equal(unname(bias$intercept), unname(coefficients[1, 1]))
  expect_equal(unname(bias$se.intercept), unname(coefficients[1, 2]))
})

test_that("longarm binary, continuous and count branches match to.long", {
  studlab <- c("Study Z", "Study A", "Study M")
  treat1 <- c("Drug Z", "Drug A", "Drug M")
  treat2 <- c("Control Z", "Control A", "Control M")
  n1 <- c(45, 62, 78)
  n2 <- c(51, 58, 83)
  event1 <- c(8, 15, 21)
  event2 <- c(12, 11, 18)
  mean1 <- c(4.2, 6.5, 8.1)
  mean2 <- c(3.7, 7.2, 6.8)
  sd1 <- c(1.1, 2.3, 1.7)
  sd2 <- c(1.8, 1.5, 2.6)
  time1 <- c(81, 124, 156)
  time2 <- c(102, 116, 166)

  binary <- longarm(treat1, treat2, studlab = studlab,
                    event1 = event1, n1 = n1, event2 = event2, n2 = n2,
                    keep.internal = TRUE)
  continuous <- longarm(treat1, treat2, studlab = studlab,
                        n1 = n1, mean1 = mean1, sd1 = sd1,
                        n2 = n2, mean2 = mean2, sd2 = sd2,
                        keep.internal = TRUE)
  count <- longarm(treat1, treat2, studlab = studlab,
                   event1 = event1, time1 = time1, n1 = n1,
                   event2 = event2, time2 = time2, n2 = n2,
                   keep.internal = TRUE)
  ref.binary <- metafor::to.long(
    "RD", ai = event1, bi = n1 - event1,
    ci = event2, di = n2 - event2, slab = studlab,
    var.names = c(".id.m4", ".grp.m4", "events", "nonevents"))
  ref.continuous <- metafor::to.long(
    "MD", n1i = n1, m1i = mean1, sd1i = sd1,
    n2i = n2, m2i = mean2, sd2i = sd2, slab = studlab,
    var.names = c(".id.m4", ".grp.m4", "mean", "sd", "n"))
  ref.count <- metafor::to.long(
    "IRD", x1i = event1, t1i = time1, x2i = event2, t2i = time2,
    slab = studlab,
    var.names = c(".id.m4", ".grp.m4", "events", "time"))

  actual <- list(binary, continuous, count)
  references <- list(ref.binary, ref.continuous, ref.count)
  arm.columns <- list(c("events", "nonevents"), c("mean", "sd", "n"),
                      c("events", "time"))
  for (i in seq_along(actual)) {
    result <- actual[[i]]
    ref <- references[[i]]
    study <- match(as.character(ref$.id.m4), studlab)
    expect_equal(nrow(result), 2L * length(studlab))
    expect_equal(result$.id.m4, ref$.id.m4)
    expect_equal(result$.grp.m4, ref$.grp.m4)
    expect_equal(as.character(result$studlab),
                 as.character(ref$.id.m4))
    expect_equal(result$treat,
                 ifelse(ref$.grp.m4 == 1,
                        treat1[study], treat2[study]))
    for (column in arm.columns[[i]])
      expect_equal(result[[column]], ref[[column]])
    expect_equal(result$n,
                 ifelse(ref$.grp.m4 == 1,
                        n1[study], n2[study]))
  }
  expect_equal(binary$n, ref.binary$events + ref.binary$nonevents)
})
