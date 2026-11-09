metabind_subgroup_fixture <- function() {
  metabin(1:9, 101:109, 2:10, 201:209,
           subgroup = rep(c("A", "B", "C"), each = 3),
           method.random.ci = rep("HK", 3),
           adhoc.hakn.ci = c("", "se", "ci"),
           prediction = TRUE,
           method.predict = c("V", "HK", "HK", "HK-PR"),
           adhoc.hakn.pi = c("", "", "se", "se"))
}

test_that("subgroup random interval keys follow the source matrix columns", {
  ma <- metabind_subgroup_fixture()
  mb <- metabind(ma, ma)
  n.by <- length(ma$subgroup.levels)
  stored <- mb$subgroup$random[mb$subgroup$random$id == 1, ]

  expect_identical(stored$method.interval,
                    rep(unname(ma$method.random.ci), each = n.by))
  expect_identical(stored$adhoc.hakn,
                    rep(unname(ma$adhoc.hakn.ci), each = n.by))
  expect_identical(stored$subgroup.levels,
                    rep(ma$subgroup.levels, length(ma$method.random.ci)))
  expect_equal(as.numeric(stored$lower), as.numeric(ma$lower.random.w))
  expect_equal(stored$Q.b, rep(unname(ma$Q.b.random), each = n.by))

  dat <- as.data.frame(mb, common = FALSE, prediction = FALSE,
                       drop_all_NA = FALSE)
  expect_equal(nrow(dat), 2 * n.by * length(ma$method.random.ci))
  expect_true(all(dat$type == "Random effects model"))
  expect_false(anyNA(dat$TE))
  first <- dat[dat$id == 1, ]
  expect_equal(as.numeric(first$lower.w), as.numeric(ma$lower.random.w))
  expect_equal(first$lower, rep(unname(ma$lower.random), each = n.by))
})

test_that("subgroup prediction keys join to the corresponding pooled intervals", {
  ma <- metabind_subgroup_fixture()
  mb <- metabind(ma, ma)
  n.by <- length(ma$subgroup.levels)
  stored <- mb$subgroup$predict[mb$subgroup$predict$id == 1, ]
  expect_identical(stored$method.interval,
                    rep(unname(ma$method.predict), each = n.by))
  expect_identical(stored$adhoc.hakn,
                    rep(unname(ma$adhoc.hakn.pi), each = n.by))
  expect_equal(as.numeric(stored$lower), as.numeric(ma$lower.predict.w))

  dat <- as.data.frame(mb, common = FALSE, random = FALSE,
                       prediction = TRUE, drop_all_NA = FALSE)
  expect_equal(nrow(dat), 2 * n.by * length(ma$method.predict))
  expect_true(all(dat$type == "Prediction interval"))
  expect_false(anyNA(dat$lower))
  first <- dat[dat$id == 1, ]
  expect_equal(as.numeric(first$lower.w), as.numeric(ma$lower.predict.w))
  expect_equal(first$lower, rep(unname(ma$lower.predict), each = n.by))
})

test_that("subgroup conversion retains joined rows and honors model selection", {
  ma <- metabind_subgroup_fixture()
  mb <- metabind(ma, ma)
  expected <- sum(vapply(mb$subgroup, nrow, integer(1)))
  dat <- as.data.frame(mb)
  expect_equal(nrow(dat), expected)
  expect_false(anyNA(dat$type))
  expect_true("subgroup.levels" %in% names(dat))

  common <- as.data.frame(mb, random = FALSE, prediction = FALSE)
  expect_equal(nrow(common), nrow(mb$subgroup$common))
  expect_true(all(common$type == "Common effect model"))

  overall <- as.data.frame(mb, subgroup = FALSE, common = FALSE,
                           prediction = FALSE)
  expect_equal(nrow(overall), nrow(mb$random))
  expect_true(all(overall$type == "Random effects model"))
  expect_false("subgroup.levels" %in% names(overall))
})

test_that("conversion without source subgroups honors model selection", {
  ma <- metagen(log(2:4), c(0.1, 0.2, 0.3), sm = "RR", prediction = TRUE)
  mb <- metabind(ma, ma)
  random <- as.data.frame(mb, common = FALSE, prediction = FALSE)
  expect_equal(nrow(random), nrow(mb$random))
  expect_true(all(random$type == "Random effects model"))
  prediction <- as.data.frame(mb, common = FALSE, random = FALSE)
  expect_equal(nrow(prediction), nrow(mb$predict))
  expect_true(all(prediction$type == "Prediction interval"))
})
