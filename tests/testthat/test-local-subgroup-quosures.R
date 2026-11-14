subgroup_quosure_data <- function() {
  data.frame(n = rep(100, 4), event = c(2, 4, 6, 8),
             event.c = c(3, 5, 7, 9), time = c(100, 200, 300, 400),
             time.c = c(110, 210, 310, 410), mean = 1:4, sd = 1,
             cor = c(0.1, 0.2, 0.3, 0.4),
             region = c("N", "N", "S", "S"),
             other = c("X", "Y", "X", "Y"))
}

test_that("constructors retain subgroup expressions and legacy precedence", {
  dat <- subgroup_quosure_data()
  inputs <- list(
    metacor = list(cor = quote(cor), n = quote(n)),
    metainc = list(event.e = quote(event), time.e = quote(time),
                   event.c = quote(event.c), time.c = quote(time.c)),
    metamean = list(n = quote(n), mean = quote(mean), sd = quote(sd)),
    metaprop = list(event = quote(event), n = quote(n), method = "Inverse"),
    metarate = list(event = quote(event), time = quote(time))
  )

  for (nm in names(inputs)) {
    fun <- get(nm, envir = asNamespace("meta"))
    args <- c(inputs[[nm]], list(data = dat))
    ma <- do.call(fun, c(args, list(subgroup = quote(region))))
    expect_identical(ma$subgroup, dat$region, info = nm)
    expect_identical(ma$subgroup.name, "region", info = nm)

    warnings <- character()
    old <- withCallingHandlers(
      do.call(fun, c(args, list(byvar = quote(region)))),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_length(warnings, 1)
    expect_identical(old$subgroup, dat$region, info = nm)
    expect_identical(old$subgroup.name, "region", info = nm)

    expect_warning(
      both <- do.call(fun, c(args, list(subgroup = quote(region),
                                       byvar = quote(other)))),
      "ignored", fixed = TRUE
    )
    expect_identical(both$subgroup, dat$region, info = nm)
    expect_identical(both$subgroup.name, "region", info = nm)

    count <- 0L
    ma <- do.call(fun, c(args, list(
      subgroup = quote({ count <<- count + 1L; region })
    )))
    expect_equal(count, 1L, info = nm)
    expect_identical(ma$subgroup, dat$region, info = nm)
  }
})

test_that("metainc preserves automatic and user-supplied pairwise groups", {
  dat <- subgroup_quosure_data()
  pw <- pairwise(
    treat = list(rep("A", 4), c("B", "B", "C", "C")),
    event = list(dat$event, dat$event.c),
    time = list(dat$time, dat$time.c), n = list(dat$n, dat$n),
    studlab = 1:4, sm = "IRR"
  )
  pw$region <- dat$region
  expect_length(unique(metainc(pw)$subgroup), 2)
  ma <- metainc(pw, subgroup = region)
  expect_identical(ma$subgroup, pw$region)
  expect_identical(ma$subgroup.name, "region")
  expect_warning(ma <- metainc(pw, byvar = region), "deprecated", fixed = TRUE)
  expect_identical(ma$subgroup, pw$region)
  expect_identical(ma$subgroup.name, "region")

  pw$region[] <- "one group"
  expect_warning(ma <- metainc(pw, byvar = region, subset = 1:2),
                 "deprecated", fixed = TRUE)
  expect_identical(ma$subgroup, rep("one group", 2))
})

test_that("update retains stored groups and captures replacement expressions", {
  dat <- subgroup_quosure_data()
  ma <- metagen(mean, sd, data = dat)
  grouped <- update(ma, subgroup = region)
  expect_identical(grouped$subgroup, dat$region)
  expect_identical(grouped$subgroup.name, "region")
  expect_identical(update(grouped, method.tau = "DL")$subgroup, dat$region)
  expect_null(update(grouped, subgroup = NULL)$subgroup)

  expect_warning(legacy <- update(ma, byvar = region), "deprecated", fixed = TRUE)
  expect_identical(legacy$subgroup, dat$region)
  expect_identical(legacy$subgroup.name, "region")
  expect_warning(both <- update(ma, subgroup = region, byvar = other),
                 "ignored", fixed = TRUE)
  expect_identical(both$subgroup, dat$region)
  expect_identical(both$subgroup.name, "region")

  count <- 0L
  grouped <- update(ma, subgroup = { count <<- count + 1L; region })
  expect_equal(count, 1L)
  expect_identical(grouped$subgroup, dat$region)
})

test_that("metacont retains subgroup expressions and evaluates them once", {
  dat <- data.frame(n.e = rep(100, 4), mean.e = 1:4, sd.e = 1,
                    n.c = rep(100, 4), mean.c = 2:5, sd.c = 1,
                    region = c("N", "N", "S", "S"))
  ma <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
                 data = dat, subgroup = region)
  expect_identical(ma$subgroup, dat$region)
  expect_identical(ma$subgroup.name, "region")

  count <- 0L
  ma <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data = dat,
                 subgroup = { count <<- count + 1L; region })
  expect_equal(count, 1L)
  expect_identical(ma$subgroup, dat$region)
})

test_that("metacont preserves pairwise defaults and explicit legacy groups", {
  pw <- pairwise(
    treat = list(rep("A", 4), c("B", "B", "C", "C")),
    n = list(rep(100, 4), rep(100, 4)),
    mean = list(1:4, 2:5), sd = list(rep(1, 4), rep(1, 4)),
    studlab = 1:4, sm = "MD"
  )
  pw$region <- c("N", "N", "S", "S")
  expect_length(unique(metacont(pw)$subgroup), 2)

  ma <- metacont(pw, subgroup = region)
  expect_identical(ma$subgroup, pw$region)
  expect_identical(ma$subgroup.name, "region")

  warnings <- character()
  ma <- withCallingHandlers(
    metacont(pw, byvar = region),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1)
  expect_match(warnings, "deprecated", fixed = TRUE)
  expect_identical(ma$subgroup, pw$region)
  expect_identical(ma$subgroup.name, "region")

  pw$region[] <- "one group"
  expect_warning(ma <- metacont(pw, byvar = region, subset = 1:2),
                 "deprecated", fixed = TRUE)
  expect_identical(ma$subgroup, rep("one group", 2))
})
