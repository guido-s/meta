test_that("method ordering and slot renumbering do not change sharing", {
  x <- interval_slot_fixture(TRUE)
  for (family in c("random", "predict")) {
    rows <- x[[family]]$id == 2L
    x[[family]]$interval.slot[rows] <- rev(x[[family]]$interval.slot[rows]) + 10L
    rows <- x$subgroup[[family]]$id == 2L
    x$subgroup[[family]]$interval.slot[rows] <-
      13L - x$subgroup[[family]]$interval.slot[rows]
  }
  expect_true(meta:::prepare_metabind_subgroups(x)$shared.overall)
  for (family in c("random", "predict")) {
    selected <- x
    selected$common <- NULL
    other <- if (family == "random") "predict" else "random"
    selected[[other]] <- NULL
    for (part in c("pooled", "subgroup")) {
      dat <- if (part == "pooled") selected[[family]] else selected$subgroup[[family]]
      dat <- dat[dat$interval.slot %in% c(2L, 11L), ]
      if (part == "pooled") selected[[family]] <- dat else
        selected$subgroup[[family]] <- dat
    }
    expect_true(meta:::prepare_metabind_subgroups(selected)$shared.overall)
    expect_silent(drawn <- interval_slot_draw(selected, backtransf = FALSE))
    expect_length(drawn$args$x[[paste0("lower.", family)]], 1L)
    expect_silent(capture.output(print(selected)))
  }
})

test_that("occurrence ranks handle missing IDs and preserve within-family multiplicity", {
  dat <- data.frame(id = c(NA_integer_, NA_integer_, 2L, 2L),
                    interval.slot = c(2L, 1L, 12L, 11L))
  descriptors <- data.frame(method = rep("HK", 4))
  expect_identical(meta:::metabind_interval_occurrence(dat, descriptors),
                   c(2L, 1L, 2L, 1L))
})
