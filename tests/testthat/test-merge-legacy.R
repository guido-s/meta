test_that("merge_legacy preserves the inferred prefix for NULL input", {
  print <- NULL
  res <- meta:::merge_legacy(
    print, list(print.tau2.ci = TRUE), list(tau2.ci = FALSE)
  )

  expect_true(res$res$tau2.ci)
  expect_length(res$args, 0)
  expect_false(res$miss$tau2.ci)
})

test_that("merge_legacy maps and removes legacy list arguments", {
  dots <- list(fs.main = 3, graphical = TRUE)

  expect_warning(
    res <- meta:::merge_legacy(list(heading = 2), dots,
                               list(main = 1, heading = 4), "fs", warn = TRUE),
    "deprecated"
  )

  expect_equal(res$res[c("main", "heading")], list(main = 3, heading = 2))
  expect_equal(res$args, list(graphical = TRUE))
  expect_false(res$miss$main)
  expect_false(res$miss$heading)
})

test_that("merge_legacy supports different prefixes", {
  expect_warning(
    ff <- meta:::merge_legacy(NULL, list(ff.study = "italic"),
                              list(study = "plain"), "ff", warn = TRUE),
    "deprecated"
  )
  expect_warning(
    col <- meta:::merge_legacy(NULL, list(col.study = "red"),
                               list(study = "black"), "col", warn = TRUE),
    "deprecated"
  )

  expect_equal(ff$res$study, "italic")
  expect_equal(col$res$study, "red")
  expect_false(ff$miss$study)
  expect_false(col$miss$study)
})

test_that("merge_legacy rejects explicitly supplied list elements", {
  expect_error(
    meta:::merge_legacy(list(main = 1), list(fs.main = 2),
                        list(main = 3), "fs"),
    "cannot both be specified"
  )
})

test_that("merge_legacy inherits unspecified child elements", {
  res <- meta:::merge_legacy(
    list(parent = "blue"),
    list(),
    list(parent = "black", child = "grey", grandchild = "white"),
    "test",
    inherit = list(child = "parent", grandchild = "child")
  )

  expect_equal(res$res[c("parent", "child", "grandchild")],
               list(parent = "blue", child = "blue", grandchild = "blue"))
  expect_true(res$miss$child)
  expect_true(res$miss$grandchild)
})

test_that("merge_legacy preserves modern and legacy named NULL values", {
  defaults <- list(parent = "blue", child = "grey")
  for (legacy in c(FALSE, TRUE)) {
    if (legacy)
      res <- meta:::merge_legacy(NULL, list(test.child = NULL), defaults, "test")
    else
      res <- meta:::merge_legacy(list(child = NULL), list(), defaults, "test")
    expect_named(res$res, names(defaults))
    expect_null(res$res$child)
    expect_false(res$miss$child)
    expect_true(res$null$child)
    expect_identical(names(res$miss), names(res$null))
  }
  res <- meta:::merge_legacy(list(child = NULL), list(), defaults, "test",
                             inherit = list(child = "parent"))
  expect_identical(res$res$child, "blue")
  expect_false(res$null$child)
})

test_that("NULL inheritance follows updated parents independently of mapping order", {
  defaults <- list(child = 1, parent = 2, grandparent = 1)
  mappings <- list(list(child = "parent", parent = "grandparent"),
                    list(parent = "grandparent", child = "parent"))
  for (mapping in mappings) {
    for (legacy in c(FALSE, TRUE)) {
      if (legacy)
        res <- meta:::merge_legacy(list(grandparent = 3), list(test.child = NULL),
                                    defaults, "test", inherit = mapping)
      else
        res <- meta:::merge_legacy(list(child = NULL, grandparent = 3), list(),
                                    defaults, "test", inherit = mapping)
      expect_identical(res$res, list(child = 3, parent = 3, grandparent = 3))
      expect_false(res$miss$child)
      expect_false(res$null$child)
      expect_identical(names(res$miss), names(res$null))
    }
  }
})
