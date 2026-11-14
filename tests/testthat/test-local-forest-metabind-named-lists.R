metabind_named_source_fixture <- function(method.tau = "DL", nonshared = FALSE) {
  dat <- data.frame(
    event.e = c(6, 30, 12, 55, 8, 45, 18, 60),
    n.e = c(100, 110, 90, 120, 100, 110, 90, 120),
    event.c = c(20, 12, 25, 18, 24, 15, 30, 20),
    n.c = c(100, 110, 90, 120, 100, 110, 90, 120),
    first = rep(c("A", "B"), each = 4),
    second = rep(c("C", "D"), 4))
  ma <- metabin(event.e, n.e, event.c, n.c, data = dat,
                subgroup = first, sm = "RR", method = "Inverse",
                method.tau = "DL", method.random.ci = "classic",
                common = TRUE, random = TRUE, prediction = FALSE)
  x <- metabind(ma, update(ma, subgroup = second, method.tau = method.tau),
                name = c("First grouping", "Second grouping"),
                keepdata = TRUE, keepmeta = TRUE)
  if (nonshared)
    x$samedata <- FALSE
  x
}

metabind_named_plain_fixture <- function() {
  ma <- metagen(c(-0.4, 0.8, 0.5, -0.3, 0.9, 0.2),
                 c(0.1, 0.3, 0.2, 0.2, 0.4, 0.15),
                 sm = "MD", method.tau = "DL",
                 common = TRUE, random = TRUE, prediction = FALSE)
  metabind(ma, update(ma, method.tau = "REML"))
}

metabind_named_draw <- function(x, ...) {
  observed <- labels <- columns <- NULL
  real.forest <- meta:::forest.meta
  real.internal <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest.meta = function(...) {
      observed <<- list(...)
      do.call(real.forest, observed)
    },
    forest_meta_internal = function(...) {
      args <- list(...)
      columns <<- args$cols
      labels <<- lapply(columns$col.studlab$labels, `[[`, "label")
      do.call(real.internal, args)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- meta::forest(x, layout = "meta", ...)
  expect_false(anyDuplicated(names(observed)) > 0)
  expect_type(observed$test, "list")
  expect_type(observed$text, "list")
  expect_false(any(grepl("^test[.]|^text[.]", names(observed))))
  if (is.list(observed$label))
    expect_false(any(grepl("^label[.]", names(observed))))
  list(result = result, args = observed, labels = labels, columns = columns)
}

metabind_named_effect_labels <- function(drawn, prefix = "Test for overall effect") {
  labels <- drawn$labels[!is.na(drawn$columns$col.studlab$rows)]
  Filter(function(z) is.call(z) &&
            grepl(prefix, paste(deparse(z[[2]]), collapse = " "), fixed = TRUE),
         labels)
}

metabind_named_label_text <- function(drawn) {
  vapply(drawn$labels, function(z) paste(deparse(z), collapse = " "), "")
}

test_that("metabind named-list formals follow forest.meta argument order", {
  nms <- names(formals(meta:::forest.metabind))
  lists <- c("test", "label", "col", "print", "type", "text", "digits")
  native <- names(formals(meta:::forest.meta))
  expect_true(all(lists %in% nms))
  expect_identical(nms[nms %in% lists], native[native %in% lists])
  expect_identical(nms[match("leftcols", nms) + 0:3],
                   c("leftcols", "rightcols", "leftlabs", "rightlabs"))
  expect_equal(match("type", nms) + 1L, match("text", nms))
  expect_equal(match("text", nms) + 1L, match("digits", nms))
})

test_that("overall named control disables shared and component effect tests", {
  for (kind in c("shared", "tau", "data")) {
    x <- metabind_named_source_fixture(
      method.tau = if (kind == "tau") "REML" else "DL",
      nonshared = kind == "data")
    baseline <- metabind_named_draw(x)
    expect_length(metabind_named_effect_labels(baseline),
                  if (kind == "shared") 2 else 4)
    expect_identical(baseline$args$overall, kind == "shared")
    expect_true(baseline$args$details)
    expect_match(paste(baseline$args$x$.text.details.methods, collapse = " "),
                 "Null hypothesis", fixed = TRUE)
    expect_warning(disabled <- metabind_named_draw(x, test = list(overall = FALSE)),
                   NA)
    expect_length(metabind_named_effect_labels(disabled), 0)
    for (field in c("overall.common", "overall.random",
                    "effect.subgroup.common", "effect.subgroup.random"))
      expect_false(any(disabled$args$test[[field]]))
    expect_false(grepl("Null hypothesis",
                       paste(disabled$args$x$.text.details.methods, collapse = " "),
                       fixed = TRUE))
    expect_true(disabled$args$details)
  }
})

test_that("model-specific overall controls match their legacy spelling", {
  for (nonshared in c(FALSE, TRUE)) {
    x <- metabind_named_source_fixture(nonshared = nonshared)
    modern <- metabind_named_draw(x, test = list(overall.common = FALSE))
    legacy <- metabind_named_draw(x, test.overall.common = FALSE)
    expect_equal(modern$args$test, legacy$args$test)
    expect_equal(modern$labels, legacy$labels)
    rows <- metabind_named_effect_labels(modern)
    expect_length(rows, if (nonshared) 2 else 1)
    expect_true(all(vapply(rows, function(z)
      grepl("random", z[[2]], fixed = TRUE), logical(1))))
  }
})

test_that("component test controls and per-component model masks are respected", {
  x <- metabind_named_source_fixture(method.tau = "REML")
  modern <- metabind_named_draw(x, test = list(effect.subgroup = FALSE))
  legacy <- metabind_named_draw(x, test.effect.subgroup = FALSE)
  expect_equal(modern$args$test, legacy$args$test)
  expect_equal(modern$labels, legacy$labels)
  expect_length(metabind_named_effect_labels(modern), 0)
  masked <- metabind_named_draw(
    x, test = list(effect.subgroup.common = c(TRUE, FALSE),
                   effect.subgroup.random = c(FALSE, TRUE)))
  expect_equal(masked$args$test$effect.subgroup.common, c(TRUE, FALSE))
  expect_equal(masked$args$test$effect.subgroup.random, c(FALSE, TRUE))
  expect_length(masked$args$x$subgroup.levels, 2)
  rows <- metabind_named_effect_labels(masked)
  expect_length(rows, 2)
  expect_match(rows[[1]][[2]], "common", fixed = TRUE)
  expect_match(rows[[2]][[2]], "random", fixed = TRUE)
  for (field in c("effect.subgroup.common", "effect.subgroup.random"))
    expect_error(metabind_named_draw(
      x, test = setNames(list(c(TRUE, FALSE, TRUE)), field)),
      "length|Length")
})

test_that("display sorting preserves the component identity of test masks", {
  components <- lapply(seq_len(11), function(i)
    metagen(i + seq_len(6) / 10, rep(0.1, 6),
             subgroup = rep(c("A", "B"), each = 3), method.tau = "DL"))
  x <- metabind(components, name = paste0("Component", seq_len(11)))
  enabled <- seq_len(11) == 2
  drawn <- metabind_named_draw(x, sort.subgroup = TRUE, details = FALSE,
                               test = list(effect.subgroup.common = enabled,
                                           effect.subgroup.random = FALSE))
  expect_identical(as.character(drawn$args$x$subgroup.levels), sort(as.character(seq_len(11))))
  expect_equal(which(drawn$args$test$effect.subgroup.common),
               match("2", as.character(drawn$args$x$subgroup.levels)))
  expect_length(metabind_named_effect_labels(drawn), 1)
})

test_that("named labels and model text render together without duplicate arguments", {
  x <- metabind_named_source_fixture(nonshared = TRUE)
  expect_warning(drawn <- metabind_named_draw(
    x, label = list(left = "Favours custom", pval = "Custom P",
                     test.effect.subgroup.common = "MY TEST"),
    text = list(common = "Custom common", random = "Custom random"),
    leftcols = c("studlab", "effect", "ci", "pval")), NA)
  expect_type(drawn$args$label, "list")
  expect_identical(drawn$args$label$left, "Favours custom")
  expect_identical(drawn$args$label$pval, "Custom P")
  expect_true("Custom P" %in% drawn$result$leftlabs)
  expect_length(metabind_named_effect_labels(drawn, "MY TEST"), 2)
  rendered <- metabind_named_label_text(drawn)
  expect_true(any(grepl("Custom common", rendered, fixed = TRUE)))
  expect_true(any(grepl("Custom random", rendered, fixed = TRUE)))
  expect_identical(drawn$args$text$common.w, "Custom common")
  expect_identical(drawn$args$text$random.w, "Custom random")
})

test_that("explicit child model text takes precedence over inherited parent text", {
  drawn <- metabind_named_draw(
    metabind_named_source_fixture(nonshared = TRUE),
    text = list(common = "Parent common", random = "Parent random",
                common.w = "Child common", random.w = "Child random"))
  expect_identical(drawn$args$text$common, "Parent common")
  expect_identical(drawn$args$text$random, "Parent random")
  expect_identical(drawn$args$text$common.w, "Child common")
  expect_identical(drawn$args$text$random.w, "Child random")
  rendered <- metabind_named_label_text(drawn)
  expect_true(any(grepl("Child common", rendered, fixed = TRUE)))
  expect_true(any(grepl("Child random", rendered, fixed = TRUE)))
})

test_that("legacy labels and text normalize once and conflicting elements fail early", {
  x <- metabind_named_source_fixture(nonshared = TRUE)
  expect_warning(drawn <- metabind_named_draw(
    x, label.e = "Custom experimental", text.common = "Legacy common",
    text.random.w = "Legacy random child",
    label.test.effect.subgroup.common = "LEGACY TEST"), NA)
  expect_identical(drawn$args$label$e, "Custom experimental")
  expect_identical(drawn$args$text$common, "Legacy common")
  expect_identical(drawn$args$text$common.w, "Legacy common")
  expect_identical(drawn$args$text$random.w, "Legacy random child")
  expect_length(metabind_named_effect_labels(drawn, "LEGACY TEST"), 2)
  conflicts <- list(
    list(test = list(overall = FALSE), test.overall = FALSE),
    list(label = list(e = "A"), label.e = "A"),
    list(text = list(common = "A"), text.common = "A"),
    list(text = list(random.w = "A"), text.random.w = "A"),
    list(label = list(test.effect.subgroup.common = "A"),
         label.test.effect.subgroup.common = "A"))
  for (args in conflicts)
    expect_error(do.call(metabind_named_draw, c(list(x = x), args)),
                 "cannot both be specified", fixed = TRUE)
})

test_that("axis label scalars and deprecated labels remain supported", {
  x <- metabind_named_source_fixture()
  expect_warning(hidden <- metabind_named_draw(x, label = FALSE), NA)
  expect_identical(hidden$args$label, FALSE)
  expect_length(metabind_named_effect_labels(hidden), 2)
  expect_warning(old <- metabind_named_draw(
    x, labels = FALSE, warn.deprecated = TRUE), "deprecated")
  expect_identical(old$args$label, FALSE)
  expect_false("labels" %in% names(old$args))
  expect_equal(old$labels, hidden$labels)
  ticks <- c("Low", "Null", "High")
  expect_warning(custom <- metabind_named_draw(
    x, label = ticks, at = c(-2, 0, 2), backtransf = FALSE), NA)
  expect_identical(custom$args$label, ticks)
  expect_equal(custom$args$at, c(-2, 0, 2))
})

test_that("no-source plots filter synthetic tests without pooling dummy results", {
  x <- metabind_named_plain_fixture()
  expect_null(x$subgroup)
  expect_warning(baseline <- metabind_named_draw(x), NA)
  expect_warning(explicit <- metabind_named_draw(
    x, test = list(overall = TRUE), label = list(left = "Custom axis"),
    text = list(common = "Custom common", random = "Custom random")), NA)
  for (drawn in list(baseline, explicit)) {
    expect_false(any(unlist(drawn$args$test)))
    expect_false(any(c("subgroup", "subgroup.common", "subgroup.random",
                       "effect.subgroup", "effect.subgroup.common",
                       "effect.subgroup.random") %in% names(drawn$args$test)))
    expect_false(drawn$args$x$overall)
    expect_false(drawn$args$x$common)
    expect_false(drawn$args$x$random)
    expect_length(metabind_named_effect_labels(drawn), 0)
    expect_false(drawn$args$details)
    expect_equal(nrow(drawn$result$data), nrow(as.data.frame(x)))
  }
  expect_equal(explicit$result$data$TE.orig, baseline$result$data$TE.orig)
  expect_identical(explicit$args$label$left, "Custom axis")
  expect_identical(explicit$args$text$common, "Custom common")
  expect_identical(explicit$args$text$random, "Custom random")
})

test_that("between-component named test requests warn and remain disabled", {
  expect_warning(drawn <- metabind_named_draw(
    metabind_named_source_fixture(), test = list(subgroup = TRUE)),
    "Between-component tests are not available", fixed = TRUE)
  for (field in c("subgroup", "subgroup.common", "subgroup.random"))
    expect_false(drawn$args$test[[field]])
  expect_length(metabind_named_effect_labels(drawn), 2)
})

test_that("overall FALSE suppresses shared tests and null details, not native details", {
  x <- metabind_named_source_fixture()
  drawn <- metabind_named_draw(x, overall = FALSE)
  expect_false(drawn$args$overall)
  expect_length(metabind_named_effect_labels(drawn), 0)
  expect_false(drawn$args$test$overall.common)
  expect_false(drawn$args$test$overall.random)
  expect_true(drawn$args$details)
  expect_false(grepl("Null hypothesis",
                     paste(drawn$args$x$.text.details.methods, collapse = " "),
                     fixed = TRUE))
  native <- metabind_named_draw(x)
  expect_true(native$args$details)
  expect_true(length(native$args$x$.text.details.methods) > 0)
  disabled <- metabind_named_draw(x, details = FALSE)
  expect_false(disabled$args$details)
  expect_null(disabled$args$x$.text.details.methods)
})

test_that("undefined proportion null effects do not produce component test rows", {
  death <- c(3, 6, 10, 1)
  animals <- c(11, 17, 21, 6)
  grp1 <- LETTERS[c(1, 2, 2, 1)]
  grp2 <- LETTERS[rep(c(1, 2), c(2, 2))]
  ma1 <- metaprop(death, animals, sm = "PFT", subgroup = grp1, common = FALSE)
  ma2 <- metaprop(death, animals, sm = "PAS", subgroup = grp2, common = FALSE)
  drawn <- metabind_named_draw(metabind(ma1, ma2))
  expect_length(metabind_named_effect_labels(drawn), 0)
  expect_false(any(drawn$args$test$effect.subgroup.random))
  expect_false(drawn$args$test$effect.subgroup)

  shared <- metabind_named_draw(metabind(ma1, update(ma1, subgroup = grp2)))
  expect_length(metabind_named_effect_labels(shared), 0)
  expect_false(shared$args$test$overall)
  expect_false(shared$args$test$overall.random)
})

test_that("null-effect eligibility is matched to each component and model", {
  death <- c(3, 6, 10, 1)
  animals <- c(11, 17, 21, 6)
  first <- metaprop(death, animals, sm = "PFT", subgroup = LETTERS[c(1, 2, 2, 1)])
  second <- metaprop(death, animals, sm = "PAS", subgroup = LETTERS[c(1, 1, 2, 2)])
  for (defined in c(1, 2)) {
    components <- list(first, second)
    components[[defined]] <- update(components[[defined]], null.effect = 0.2)
    x <- metabind(components)
    x$subgroup$common <- x$subgroup$common[order(-x$subgroup$common$id), ]
    x$subgroup$random <- x$subgroup$random[order(-x$subgroup$random$id), ]
    for (sorted in c(FALSE, TRUE)) {
      drawn <- metabind_named_draw(x, sort.subgroup = sorted)
      expected <- as.character(drawn$args$x$subgroup.levels) == as.character(defined)
      expect_equal(drawn$args$test$effect.subgroup.common, expected)
      expect_equal(drawn$args$test$effect.subgroup.random, expected)
      expect_length(metabind_named_effect_labels(drawn), 2)
      expect_match(paste(drawn$args$x$.text.details.methods, collapse = " "),
                   "Null hypothesis", fixed = TRUE)
    }
  }
})
