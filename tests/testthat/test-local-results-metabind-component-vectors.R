test_that("variance list columns retain complete source vectors in every result row", {
  for (subgroup in c(FALSE, TRUE)) {
    sources <- metabind_component_fixture(subgroup)
    x <- metabind(sources[[1]], sources[[2]], keepmeta = TRUE)
    before <- x
    expect_length(sources[[1]]$tau2, 1L)
    expect_length(sources[[2]]$tau2, 2L)
    tables <- list(x$overall, as.data.frame(x, subgroup = FALSE, drop_all_NA = FALSE))
    if (subgroup)
      tables <- c(tables, list(as.data.frame(x, subgroup = TRUE, drop_all_NA = FALSE)))
    for (dat in tables)
      for (field in metabind_component_fields) {
        expect_type(dat[[field]], "list")
        for (row in seq_len(nrow(dat)))
          expect_identical(dat[[field]][[row]], sources[[dat$id[row]]][[field]],
                           info = field)
      }
    if (subgroup) {
      prepared <- meta:::prepare_metabind_subgroups(x)
      expect_false(prepared$shared.overall)
      for (dat in list(prepared$dat, prepared$subdat))
        for (field in metabind_component_fields)
          for (row in seq_len(nrow(dat)))
            expect_identical(dat[[field]][[row]], sources[[dat$id[row]]][[field]])
    }
    expect_silent(output <- capture.output(print(
      x, print = list(tau2 = TRUE, tau = TRUE), details.methods = FALSE)))
    expected <- meta:::format_metabind_variance(list(sources[[2]]$tau2), digits = 4)
    if (!subgroup)
      expect_true(any(grepl(expected, output, fixed = TRUE)))
    expect_identical(x, before)
    expect_identical(x$meta, sources)
  }
})

test_that("native forest totals are local and have no aggregate component CI", {
  sources <- metabind_component_fixture(TRUE)
  x <- metabind(sources[[1]], sources[[2]], keepmeta = TRUE)
  before <- x
  observed <- NULL
  draw <- meta:::forest.meta
  testthat::local_mocked_bindings(forest.meta = function(...) {
    observed <<- list(...)
    do.call(draw, observed)
  }, .package = "meta")
  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(forest(x, layout = "meta", details = FALSE,
                      print = list(tau2 = TRUE, tau2.ci = TRUE,
                                   tau = TRUE, tau.ci = TRUE)))
  native <- observed$x
  expect_equal(unname(native$tau2.w), vapply(sources, function(ma) sum(ma$tau2), 0))
  expect_equal(unname(native$tau.w),
               vapply(sources, function(ma) sqrt(sum(ma$tau^2)), 0))
  for (field in c("lower.tau2", "upper.tau2", "se.tau2", "lower.tau", "upper.tau")) {
    expect_identical(unname(native[[paste0(field, ".w")]][2]), NA_real_, info = field)
    expect_equal(unname(native[[paste0(field, ".w")]][1]),
                 unname(sources[[1]][[field]]), info = field)
  }
  expect_identical(x, before)
  expect_identical(x$meta, sources)

  sources <- metabind_component_fixture(FALSE)
  x <- metabind(sources[[1]], sources[[2]], keepmeta = TRUE)
  before <- x
  expect_silent(forest(x, layout = "meta", details = FALSE,
                      print = list(tau2 = TRUE, tau2.ci = TRUE,
                                   tau = TRUE, tau.ci = TRUE)))
  expect_identical(x, before)
  expect_identical(x$meta, sources)
})
