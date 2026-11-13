metamerge_method_fixture <- function(methods, method = "GLMM", common = FALSE) {
  dat <- data.frame(
    event.e = c(10, 20, 30, 40, 50, 30, 40, 50, 60, 70),
    event.c = c(20, 22, 25, 30, 35, 20, 25, 30, 35, 40),
    total = 100, group = rep(c("A", "B"), each = 5))
  metabin(event.e, total, event.c, total, data = dat, method = method,
             sm = "OR", common = common, subgroup = group, tau.common = TRUE,
             method.random.ci = methods)
}

test_that("unequal-method merges retain each source subgroup matrix and DF block", {
  single <- metamerge_method_fixture("HK")
  multiple <- metamerge_method_fixture(c("classic", "HK"))
  reversed <- metamerge_method_fixture(c("HK", "classic"))
  for (sources in list(list(single, multiple), list(multiple, single),
                       list(multiple, reversed))) {
    merged <- metamerge(sources[[1]], sources[[2]], label1 = "Left", label2 = "Right")
    estimate <- meta:::list2mat(merged$TE.random.w, merged$lower.random.w)
    lower <- meta:::list2mat(merged$lower.random.w)
    upper <- meta:::list2mat(merged$upper.random.w)
    methods <- unlist(lapply(sources, `[[`, "method.random.ci"), use.names = FALSE)
    expect_equal(dim(estimate), c(2L, length(methods)))
    expect_equal(dim(lower), dim(estimate))
    expect_identical(dimnames(estimate), dimnames(lower))
    expect_identical(dimnames(upper), dimnames(lower))
    expect_identical(unname(merged$method.random.ci), methods)
    expected.df <- list()
    column <- 0L
    for (source in seq_along(sources)) {
      ma <- sources[[source]]
      for (method in seq_along(ma$method.random.ci)) {
        column <- column + 1L
        for (g in seq_along(ma$subgroup.levels)) {
          point <- if (is.matrix(ma$TE.random.w)) ma$TE.random.w[g, method] else ma$TE.random.w[g]
          low <- if (is.matrix(ma$lower.random.w)) ma$lower.random.w[g, method] else ma$lower.random.w[g]
          upp <- if (is.matrix(ma$upper.random.w)) ma$upper.random.w[g, method] else ma$upper.random.w[g]
          expect_equal(unname(estimate[g, column]), unname(point))
          expect_equal(unname(lower[g, column]), unname(low))
          expect_equal(unname(upper[g, column]), unname(upp))
        }
        expected.df <- c(expected.df,
                          if (is.list(ma$df.Q.b.random)) ma$df.Q.b.random[method]
                          else list(ma$df.Q.b.random))
      }
    }
    expect_identical(unname(merged$df.Q.b.random), unname(expected.df))
    expect_equal(unname(merged$Q.b.random),
                  unlist(lapply(sources, `[[`, "Q.b.random"), use.names = FALSE))
    expect_equal(unname(merged$pval.Q.b.random),
                  unlist(lapply(sources, `[[`, "pval.Q.b.random"), use.names = FALSE))
  }
})

expect_metamerge_overall_alignment <- function(merged, sources, model, labels,
                                               output, selected = c(TRUE, TRUE)) {
  fields <- c("TE", "seTE", "lower", "upper", "statistic", "pval",
               if (model == "random") "df")
  expected <- list()
  for (source in which(selected)) {
    ma <- sources[[source]]
    n.result <- length(ma[[paste0("lower.", model)]])
    for (j in seq_len(n.result)) {
      row <- lapply(fields, function(field) {
        value <- ma[[paste0(field, ".", model)]]
        if (is.null(value)) NA_real_ else unname(value[if (length(value) == 1L) 1L else j])
      })
      names(row) <- fields
      k <- ma$k
      offset <- if (model == "random") length(ma$lower.common) else 0L
      row$k <- unname(k[if (length(k) == 1L) 1L else offset + j])
      text <- ma[[paste0("text.", model)]]
      row$text <- paste0(text[if (length(text) == 1L) 1L else j], " (", labels[source], ")")
      expected[[length(expected) + 1L]] <- as.data.frame(row)
    }
  }
  expected <- do.call(rbind, expected)
  for (field in fields) {
    actual <- merged[[paste0(field, ".", model)]]
    expect_length(actual, nrow(expected))
    expect_equal(unname(actual), expected[[field]], info = paste(model, field))
  }
  expect_identical(unname(merged[[paste0("text.", model)]]), expected$text)
  offset <- if (model == "random") length(merged$lower.common) else 0L
  actual.k <- if (length(merged$k) == 1L) rep(merged$k, nrow(expected))
  else merged$k[offset + seq_len(nrow(expected))]
  expect_equal(unname(actual.k), expected$k)
  for (j in seq_len(nrow(expected))) {
    line <- output[startsWith(output, expected$text[j])]
    expect_length(line, 1)
    for (field in c("TE", "lower", "upper", "statistic"))
      expect_match(line, trimws(meta:::formatN(expected[[field]][j], digits = 8)),
                   fixed = TRUE)
    expect_match(line, trimws(meta:::formatPT(expected$pval[j], digits = 8,
                                              scientific = FALSE, zero = TRUE,
                                              JAMA = FALSE)), fixed = TRUE)
  }
  invisible(expected)
}

test_that("overall normalization expands shared fields and preserves per-method results", {
  ma <- metamerge_method_fixture(c("classic", "HK"))
  original <- ma
  expect_length(ma$TE.random, 1)
  normalized <- meta:::normalize_overall_results(ma)
  expect_identical(normalized$TE.random, rep(ma$TE.random, 2))
  for (field in c("seTE.random", "lower.random", "upper.random", "statistic.random",
                  "pval.random", "df.random", "text.random", "method.random.ci"))
    expect_identical(normalized[[field]], ma[[field]])
  expect_identical(normalized$df.Q.b.random, ma$df.Q.b.random)
  expect_identical(meta:::normalize_overall_results(normalized), normalized)
  expect_identical(ma, original)

  missing <- ma
  missing$pval.random <- missing$df.random <- missing$text.random <- NULL
  normalized <- meta:::normalize_overall_results(missing)
  expect_identical(normalized$pval.random, rep(NA_real_, 2))
  expect_identical(normalized$df.random, rep(NA_real_, 2))
  expect_identical(normalized$text.random, rep(NA_character_, 2))
  invalid <- ma
  invalid$upper.random <- invalid$upper.random[1]
  expect_error(meta:::normalize_overall_results(invalid), "different lengths")
  invalid <- ma
  invalid$TE.random <- 1:3
  expect_error(meta:::normalize_overall_results(invalid), "does not match")
})

test_that("default overall printing aligns unequal CI counts with distinct source estimates", {
  withr::local_options(width = 300)
  single <- metamerge_method_fixture("HK", method = "Inverse")
  multiple <- metamerge_method_fixture(c("classic", "HK"))
  left.multiple <- metamerge_method_fixture(c("classic", "HK"), method = "Inverse")
  reversed <- metamerge_method_fixture(c("HK", "classic"))
  expect_gt(abs(single$TE.random - multiple$TE.random), 1e-5)
  for (sources in list(list(single, multiple), list(multiple, single),
                       list(left.multiple, reversed))) {
    merged <- metamerge(sources[[1]], sources[[2]], label1 = "Left", label2 = "Right")
    output <- capture.output(print(
      merged, details.methods = FALSE, backtransf = FALSE,
      digits = 8, digits.stat = 8, digits.pval = 8,
      scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE))
    expect_metamerge_overall_alignment(merged, sources, "random", c("Left", "Right"), output)
    expect_equal(sum(grepl("^Test for subgroup differences", output)),
                  sum(vapply(sources, function(x) length(x$method.random.ci), integer(1))))
  }
})

for (selection in list(
  list(common1 = TRUE, common2 = TRUE, random1 = TRUE, random2 = TRUE),
  list(common1 = TRUE, common2 = TRUE, random1 = FALSE, random2 = FALSE),
  list(common1 = TRUE, common2 = FALSE, random1 = FALSE, random2 = TRUE),
  list(common1 = FALSE, common2 = TRUE, random1 = TRUE, random2 = FALSE))) {
  test_that(paste("overall normalization respects selection",
                  paste(names(selection)[unlist(selection)], collapse = ", ")), {
    withr::local_options(width = 300)
    sources <- list(metamerge_method_fixture("HK", method = "Inverse", common = TRUE),
                      metamerge_method_fixture(c("classic", "HK"), common = TRUE))
    merged <- do.call(metamerge, c(list(meta1 = sources[[1]], meta2 = sources[[2]],
                                       label1 = "Left", label2 = "Right"), selection))
    output <- capture.output(print(
      merged, details.methods = FALSE, backtransf = FALSE,
      digits = 8, digits.stat = 8, digits.pval = 8,
      scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE))
    for (model in c("common", "random")) {
      selected <- unlist(selection[paste0(model, 1:2)], use.names = FALSE)
      expect_identical(merged[[model]], any(selected))
      if (any(selected))
        expect_metamerge_overall_alignment(
          merged, sources, model, c("Left", "Right"), output, selected)
      else
        expect_false(any(startsWith(output,
                                    if (model == "common") "Common effect model"
                                    else "Random effects model")))
    }
  })
}

test_that("nested merges preserve normalized overall slots and source labels", {
  withr::local_options(width = 300)
  first <- metamerge_method_fixture("HK", method = "Inverse", common = TRUE)
  second <- metamerge_method_fixture(c("classic", "HK"), common = TRUE)
  third <- metamerge_method_fixture("classic", method = "Inverse", common = TRUE)
  combined <- metamerge(first, second, label1 = "First", label2 = "Second")
  merged <- metamerge(combined, third, label1 = "Combined", label2 = "Third")
  output <- capture.output(print(
    merged, details.methods = FALSE, backtransf = FALSE,
    digits = 8, digits.stat = 8, digits.pval = 8,
    scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE))
  for (model in c("common", "random"))
    expect_metamerge_overall_alignment(
      merged, list(combined, third), model, c("Combined", "Third"), output)
  expect_equal(sum(grepl("^Test for subgroup differences", output)), 7)
})

test_that("CI-slot counts keep different study counts aligned through nested merges", {
  withr::local_options(width = 300)
  first <- update(metamerge_method_fixture("HK", method = "Inverse", common = TRUE),
                    exclude = c(FALSE, TRUE, rep(FALSE, 8)))
  second <- metamerge_method_fixture(c("classic", "HK"), common = TRUE)
  third <- metamerge_method_fixture("classic", method = "Inverse", common = TRUE)
  expect_equal(first$k, 9)
  expect_equal(second$k, 10)
  merged <- metamerge(first, second, label1 = "First", label2 = "Second")
  expect_equal(unname(merged$k), c(9, 10, 9, 10, 10))
  nested <- metamerge(merged, third, label1 = "Combined", label2 = "Third")
  expect_equal(unname(nested$k), c(9, 10, 10, 9, 10, 10, 10))
  for (case in list(list(x = merged, sources = list(first, second), labels = c("First", "Second")),
                     list(x = nested, sources = list(merged, third), labels = c("Combined", "Third")))) {
    output <- capture.output(print(
      case$x, details.methods = FALSE, backtransf = FALSE,
      digits = 8, digits.stat = 8, digits.pval = 8,
      scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE))
    for (model in c("common", "random"))
      expect_metamerge_overall_alignment(case$x, case$sources, model, case$labels, output)
  }
})

test_that("merged subgroup printing retains all methods, labels, intervals, and tests", {
  single <- metamerge_method_fixture("HK")
  multiple <- metamerge_method_fixture(c("classic", "HK"))
  reversed <- metamerge_method_fixture(c("HK", "classic"))
  for (sources in list(list(single, multiple), list(multiple, single),
                       list(multiple, reversed))) {
    merged <- metamerge(sources[[1]], sources[[2]], label1 = "Left", label2 = "Right")
    # Overall-result normalization is a separate issue from subgroup conversion.
    output <- capture.output(print(
      merged, overall = FALSE, details.methods = FALSE, backtransf = FALSE,
      digits = 8, digits.Q = 8, digits.pval.Q = 8,
      scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE))
    tables <- grep("^Results for subgroups", output)
    tests <- grep("^Test for subgroup differences", output)
    counts <- vapply(sources, function(x) length(x$method.random.ci), integer(1))
    expect_length(tables, sum(counts))
    expect_length(tests, sum(counts))
    expect_equal(sum(grepl("^Between groups", output)), sum(counts))
    labels <- rep(c("Left", "Right"), counts)
    expect_true(all(vapply(seq_along(tables), function(i)
      grepl(labels[i], output[tables[i]], fixed = TRUE), logical(1))))
    column <- 0L
    for (ma in sources) {
      for (method in seq_along(ma$method.random.ci)) {
        column <- column + 1L
        for (g in seq_along(ma$subgroup.levels)) {
          row <- output[tables[column] + 1L + g]
          values <- c(
            if (is.matrix(ma$TE.random.w)) ma$TE.random.w[g, method] else ma$TE.random.w[g],
            if (is.matrix(ma$lower.random.w)) ma$lower.random.w[g, method] else ma$lower.random.w[g],
            if (is.matrix(ma$upper.random.w)) ma$upper.random.w[g, method] else ma$upper.random.w[g])
          for (value in values)
            expect_match(row, trimws(meta:::formatN(value, digits = 8)), fixed = TRUE)
        }
        df <- if (is.list(ma$df.Q.b.random)) ma$df.Q.b.random[[method]] else ma$df.Q.b.random
        header <- output[tests[column] + 1L]
        expect_match(header, if (length(df) == 2L) "F" else "Q", fixed = TRUE)
        expected <- paste(
          "Between groups", meta:::formatN(ma$Q.b.random[method], digits = 8),
          paste(trimws(meta:::formatN(df, digits = 0)), collapse = ", "),
          meta:::formatPT(ma$pval.Q.b.random[method], digits = 8,
                         scientific = FALSE, zero = TRUE, JAMA = FALSE))
        expect_identical(
          strsplit(trimws(output[tests[column] + 2L]), "[[:space:]]+")[[1]],
          strsplit(trimws(expected), "[[:space:]]+")[[1]])
      }
    }
    summary.output <- capture.output(print(summary(merged), overall = FALSE,
                                             details.methods = FALSE))
    expect_equal(sum(grepl("^Test for subgroup differences", summary.output)), sum(counts))
  }
})

test_that("missing common DF repeat across merged result slots even with tests hidden", {
  ma <- metamerge_method_fixture("classic", method = "Inverse", common = TRUE)
  merged <- metamerge(ma, ma)
  nested <- metamerge(merged, ma)
  for (x in list(merged, nested)) {
    expect_length(x$df.Q.b, 1L)
    expect_gt(length(x$TE.common), 1L)
    for (visible in c(TRUE, FALSE)) {
      expected <- capture.output(print(
        x, random = FALSE, details.methods = FALSE, test.subgroup = visible))
      for (missing.p in c(FALSE, TRUE)) {
        missing <- x
        missing$df.Q.b.common <- NULL
        if (missing.p)
          missing$pval.Q.b.common <- NULL
        before <- missing
        expect_silent(actual <- capture.output(print(
          missing, random = FALSE, details.methods = FALSE, test.subgroup = visible)))
        expect_identical(actual, expected)
        expect_identical(missing, before)
      }
    }
  }
})
