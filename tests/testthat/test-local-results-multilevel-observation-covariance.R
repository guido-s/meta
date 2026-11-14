test_that("joint common-tau subgroup MLM filters observations consistently", {
  skip_if_not_installed("metadat")
  data(dat.konstantopoulos2011, package = "metadat")
  dat <- dat.konstantopoulos2011
  dat$se <- sqrt(dat$vi)
  dat$label <- paste0("effect-", seq_len(nrow(dat)))
  dat$group <- factor(ifelse(dat$district <= median(dat$district), "A", "B"))
  dat$exclude <- seq_len(nrow(dat)) %in% c(2, 17)
  dat$yi[8] <- Inf
  dat$se[26] <- Inf
  dat$yi[37] <- NA_real_
  dat$se[45] <- NA_real_
  keep <- !dat$exclude & is.finite(dat$yi) & is.finite(dat$se)
  retained <- dat[keep, ]
  retained$idx <- which(keep)

  ma <- metagen(yi, se, studlab = label, data = dat,
                cluster = district, rho = 0.7, exclude = exclude,
                subgroup = group, tau.common = TRUE, common = FALSE,
                method.tau = "REML", method.random.ci = "HK",
                keeprma = TRUE)
  filtered <- metagen(yi, se, studlab = label, data = retained,
                      cluster = district, rho = 0.7,
                      subgroup = group, tau.common = TRUE, common = FALSE,
                      method.tau = "REML", method.random.ci = "HK",
                      keeprma = TRUE)

  expect_identical(which(!keep), c(2L, 8L, 17L, 26L, 37L, 45L))
  expect_true(ma$three.level)
  expect_true(ma$tau.common)
  expect_identical(ma$studlab[keep], filtered$studlab)
  expect_identical(ma$cluster[keep], filtered$cluster)
  expect_identical(ma$data$.idx[keep], retained$idx)
  expect_equal(as.numeric(ma$rma.three.level$yi), as.numeric(retained$yi))
  expect_identical(as.character(ma$rma.three.level$mf.s[[1]]),
                   as.character(retained$district))
  expect_identical(as.character(ma$rma.three.level$mf.s[[2]]),
                   paste(retained$district, retained$idx, sep = "/"))
  expect_equal(ma$rma.three.level$V, filtered$rma.three.level$V)

  fields <- c("k", "k.study", "k.w", "TE.random", "seTE.random",
              "lower.random", "upper.random", "df.random", "pval.random",
              "tau2", "subgroup.levels", "TE.random.w", "seTE.random.w",
              "lower.random.w", "upper.random.w", "statistic.random.w",
              "df.random.w", "pval.random.w", "tau2.w",
              "Q.b.random", "df.Q.b.random", "pval.Q.b.random")
  expect_equal(ma[fields], filtered[fields])

  V <- metafor::vcalc(vi = retained$se^2, cluster = retained$district,
                      obs = retained$idx, rho = 0.7)
  subgroup <- metafor::rma.mv(yi, V = V, mods = ~ group - 1,
                             random = ~ 1 | district / idx, data = retained,
                             method = "REML", test = "t",
                             level = ma$level.ma * 100)
  difference <- metafor::rma.mv(yi, V = V, mods = ~ group,
                               random = ~ 1 | district / idx, data = retained,
                               method = "REML", test = "t",
                               level = ma$level.ma * 100)
  expect_identical(names(ma$TE.random.w), levels(retained$group))
  expect_equal(unname(ma$TE.random.w), as.numeric(coef(subgroup)))
  expect_equal(unname(ma$seTE.random.w), subgroup$se)
  expect_equal(unname(ma$lower.random.w), subgroup$ci.lb)
  expect_equal(unname(ma$upper.random.w), subgroup$ci.ub)
  expect_equal(unname(ma$statistic.random.w), subgroup$zval)
  expect_equal(unname(ma$df.random.w), unname(subgroup$ddf))
  expect_equal(unname(ma$pval.random.w), subgroup$pval)
  expect_equal(unname(ma$tau2.w), rep(sum(subgroup$sigma2), 2))
  expect_equal(unname(ma$Q.b.random), difference$QM)
  expect_equal(unname(ma$df.Q.b.random), difference$QMdf)
  expect_equal(unname(ma$df.Q.b.random), c(1, sum(keep) - 2))
  expect_equal(unname(ma$pval.Q.b.random), difference$QMp)
  expect_equal(unname(ma$pval.Q.b.random),
               pf(ma$Q.b.random, 1, sum(keep) - 2, lower.tail = FALSE))
})

test_that("rho-aware multilevel metareg matches direct covariance and inference", {
  skip_if_not_installed("metadat")
  data(dat.konstantopoulos2011, package = "metadat")
  dat <- dat.konstantopoulos2011
  dat$se <- sqrt(dat$vi)
  dat$group <- factor(ifelse(dat$district <= median(dat$district), "A", "B"))
  methods <- c("HK", "CR0", "CR1")
  if (suppressPackageStartupMessages(requireNamespace("clubSandwich", quietly = TRUE)))
    methods <- c(methods, "CR2")

  for (rho in c(0, 0.7)) {
    for (excluded in list(integer(), c(2L, 17L, 37L))) {
      dat$exclude <- seq_len(nrow(dat)) %in% excluded
      retained <- dat[!dat$exclude, ]
      retained$idx <- seq_len(nrow(retained))
      ma <- metagen(yi, se, data = dat, cluster = district, rho = rho,
                    exclude = exclude, common = FALSE, method.tau = "REML")
      V <- metafor::vcalc(vi = retained$se^2, cluster = retained$district,
                          obs = retained$idx, rho = rho)
      direct <- metafor::rma.mv(yi, V = V, mods = ~ group,
                                random = ~ 1 | district / idx, data = retained,
                                method = ma$method.tau, test = "z",
                                level = ma$level.ma * 100)
      direct.hk <- metafor::rma.mv(yi, V = V, mods = ~ group,
                                   random = ~ 1 | district / idx, data = retained,
                                   method = ma$method.tau, test = "t",
                                   level = ma$level.ma * 100)

      for (method in methods) {
        mr <- metareg(ma, ~ group, method.random.ci = method)
        reference <- if (method == "HK") direct.hk else
          suppressPackageStartupMessages(
            metafor::robust(direct, cluster = retained$district,
                            adjust = method != "CR0",
                            clubSandwich = method == "CR2"))
        context <- paste("rho", rho, "excluded", length(excluded), method)
        expect_s3_class(mr, "metareg")
        expect_identical(mr$.meta$method.random.ci, method)
        expect_equal(as.numeric(coef(mr)), as.numeric(coef(reference)),
                     info = context)
        fields <- c("se", "vb", "zval", "pval", "ddf", "QM", "QMdf", "QMp",
                    "sigma2", "ci.lb", "ci.ub")
        expect_equal(mr[fields], reference[fields], info = context)
        expect_equal(unname(vcov(mr)), unname(vcov(reference)), info = context)
        expect_equal(mr$V, V, info = context)
        expect_equal(as.numeric(mr$yi), as.numeric(retained$yi), info = context)
        expect_identical(as.character(mr$mf.s[[1]]),
                         as.character(retained$district))
        expect_identical(as.character(mr$mf.s[[2]]),
                         paste(retained$district, retained$idx, sep = "/"))
      }
    }
  }
})
