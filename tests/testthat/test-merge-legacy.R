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
    inherit = c(child = "parent", grandchild = "child")
  )

  expect_equal(res$res[c("parent", "child", "grandchild")],
               list(parent = "blue", child = "blue", grandchild = "blue"))
  expect_true(res$miss$child)
  expect_true(res$miss$grandchild)
})
