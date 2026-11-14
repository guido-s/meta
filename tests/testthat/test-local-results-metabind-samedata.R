metabind_samedata_sources <- function(subgroups = TRUE) {
  dat <- data.frame(
    TE = c(-0.7, 0.2, 0.8, 1.4, -0.3, 0.5, 1.1, 1.8),
    seTE = c(0.15, 0.3, 0.2, 0.4, 0.25, 0.18, 0.35, 0.22),
    study = paste0("Study ", 1:8),
    unused = 1:8,
    first = rep(c("A", "B"), each = 4),
    second = rep(c("C", "D"), 4),
    third = c("E", "F", "F", "E", "E", "F", "F", "E"))
  ma <- metagen(TE, seTE, studlab = study, data = dat, sm = "MD",
                 method.tau = "DL", method.random.ci = "classic",
                 common = TRUE, random = TRUE, prediction = TRUE,
                 method.predict = "V", keepdata = TRUE)
  if (!subgroups)
    return(list(ma, ma, ma))
  list(update(ma, subgroup = first), update(ma, subgroup = second),
       update(ma, subgroup = third))
}

test_that("source identity and retention are independent with or without subgroups", {
  for (subgroups in c(FALSE, TRUE)) {
    sources <- metabind_samedata_sources(subgroups)
    for (same in c(FALSE, TRUE)) {
      inputs <- sources
      if (!same)
        inputs[[2]]$data$unused[1] <- 99L
      # An unused source column changes identity, not any fitted estimate.
      for (field in c("TE", "seTE", "TE.common", "seTE.common",
                      "TE.random", "seTE.random", "lower.predict"))
        expect_equal(inputs[[1]][[field]], inputs[[2]][[field]])
      for (keepdata in c(FALSE, TRUE)) for (keepmeta in c(FALSE, TRUE)) {
        x <- do.call(metabind, c(inputs, list(
          keepdata = keepdata, keepmeta = keepmeta)))
        expect_true("samedata" %in% names(x))
        expect_identical(x$samedata, same)
        expect_identical("data" %in% names(x), keepdata)
        expect_identical("meta" %in% names(x), keepmeta)
        if (keepmeta)
          expect_identical(x$meta, inputs)
        if (keepdata && !same) {
          expect_length(x$data, 3L)
          for (i in seq_along(inputs))
            expect_identical(x$data[[i]], inputs[[i]]$data)
        }
        if (keepdata && same) {
          expect_s3_class(x$data, "data.frame")
          if (subgroups) {
            expect_false(".subgroup" %in% names(x$data))
            for (i in seq_along(inputs))
              expect_identical(x$data[[paste0(".subgroup", i)]],
                               inputs[[i]]$data$.subgroup)
            retained <- x$data
            retained[paste0(".subgroup", 1:3)] <- NULL
            original <- inputs[[1]]$data
            original$.subgroup <- NULL
            expect_identical(retained, original)
          }
          else
            expect_identical(x$data, inputs[[1]]$data)
        }
        expect_identical(is.null(x$subgroup), !subgroups)
        if (subgroups) {
          prepared <- meta:::prepare_metabind_subgroups(x)
          expect_identical(prepared$same.data, same)
          expect_identical(prepared$shared.overall, same)
          expect_identical(prepared$shared.heterogeneity, same)
        }
      }
    }
  }
})

test_that("absent source datasets and placeholders never prove shared data", {
  for (subgroups in c(FALSE, TRUE)) {
    sources <- metabind_samedata_sources(subgroups)[1:2]
    for (missing in c("one", "both", "placeholders")) {
      inputs <- sources
      inputs[[1]]$data <- NULL
      if (missing == "both")
        inputs[[2]]$data <- NULL
      if (missing == "placeholders") {
        inputs[[1]]$data <- data.frame(.not_a_dataset = 1)
        inputs[[2]]$data <- data.frame(.not_a_dataset = 1)
      }
      expect_equal(inputs[[1]]$TE.common, inputs[[2]]$TE.common)
      expect_equal(inputs[[1]]$TE.random, inputs[[2]]$TE.random)
      for (keepdata in c(FALSE, TRUE)) for (keepmeta in c(FALSE, TRUE)) {
        x <- do.call(metabind, c(inputs, list(
          keepdata = keepdata, keepmeta = keepmeta)))
        expect_identical(x$samedata, FALSE)
        expect_identical("data" %in% names(x), keepdata)
        expect_identical("meta" %in% names(x), keepmeta)
        if (subgroups) {
          prepared <- meta:::prepare_metabind_subgroups(x)
          expect_false(prepared$same.data)
          expect_false(prepared$shared.overall)
          expect_false(prepared$shared.heterogeneity)
          expect_true(all(c("id", "name") %in% names(prepared$dat)))
        }
      }
    }
  }
})

test_that("shared source identity does not override selected CI or tau settings", {
  sources <- metabind_samedata_sources()
  alternatives <- list(
    CI = update(sources[[2]], method.random.ci = "HK"),
    tau = update(sources[[2]], method.tau = "REML"))
  for (difference in names(alternatives)) {
    for (keepdata in c(FALSE, TRUE)) for (keepmeta in c(FALSE, TRUE)) {
      x <- metabind(sources[[1]], alternatives[[difference]],
                     keepdata = keepdata, keepmeta = keepmeta)
      expect_identical(x$samedata, TRUE)
      prepared <- meta:::prepare_metabind_subgroups(x)
      expect_true(prepared$same.data)
      expect_false(prepared$shared.overall)
      expect_identical(prepared$shared.heterogeneity, difference == "CI")
      # Only the selected common effects agree; tau metadata remains separate.
      x$random <- x$predict <- NULL
      prepared <- meta:::prepare_metabind_subgroups(x)
      expect_true(prepared$same.data)
      expect_true(prepared$shared.overall)
      expect_identical(prepared$shared.heterogeneity, difference == "CI")
      expect_equal(nrow(prepared$dat), 1L)
    }
  }
})

test_that("print and real forest share stored summaries without retained data", {
  sources <- metabind_samedata_sources()[1:2]
  selections <- list(
    common = list(common = TRUE, random = FALSE, prediction = FALSE),
    random = list(common = FALSE, random = TRUE, prediction = FALSE),
    predict = list(common = FALSE, random = FALSE, prediction = TRUE),
    all = list(common = TRUE, random = TRUE, prediction = TRUE))
  for (keepmeta in c(FALSE, TRUE)) {
    retained <- metabind(sources[[1]], sources[[2]],
                          keepdata = TRUE, keepmeta = keepmeta)
    discarded <- metabind(sources[[1]], sources[[2]],
                           keepdata = FALSE, keepmeta = keepmeta)
    # Distinguish stored results from anything freshly pooled from source rows.
    for (family in c("common", "random", "predict")) {
      for (field in intersect(c("TE", "lower", "upper"),
                              names(retained[[family]]))) {
        retained[[family]][[field]] <- retained[[family]][[field]] + 0.137
        discarded[[family]][[field]] <- retained[[family]][[field]]
      }
    }
    for (selection in selections) {
      before <- discarded
      expect_silent(printed.retained <- capture.output(do.call(
        print, c(list(x = retained, backtransf = FALSE), selection))))
      expect_silent(printed.discarded <- capture.output(do.call(
        print, c(list(x = discarded, backtransf = FALSE), selection))))
      expect_identical(printed.discarded, printed.retained)
      expect_gt(length(printed.discarded), 0L)
      expect_silent(drawn.retained <- do.call(
        metabind_render_pdf,
        c(list(x = retained, backtransf = FALSE), selection)))
      expect_silent(drawn.discarded <- do.call(
        metabind_render_pdf,
        c(list(x = discarded, backtransf = FALSE), selection)))
      expect_true(drawn.discarded$args$overall)
      expect_identical(drawn.discarded$result$data, drawn.retained$result$data)
      expect_identical(drawn.discarded$labels, drawn.retained$labels)
      for (flag in c("common", "random", "prediction", "overall",
                     "overall.hetstat", "subgroup.hetstat"))
        expect_identical(drawn.discarded$args[[flag]],
                         drawn.retained$args[[flag]])
      ma <- drawn.discarded$args$x
      for (family in c("common", "random", "predict")) {
        selected <- selection[[if (family == "predict") "prediction" else family]]
        if (!selected)
          next
        for (field in intersect(c("TE", "seTE", "lower", "upper", "df"),
                                names(discarded[[family]]))) {
          name <- paste0(field, ".", family)
          expect_equal(as.numeric(ma[[name]]), as.numeric(discarded[[family]][[field]][1]))
          expect_identical(ma[[name]], drawn.retained$args$x[[name]])
        }
      }
      expect_identical(discarded, before)
    }
  }
})
