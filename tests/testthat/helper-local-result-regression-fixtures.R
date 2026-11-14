metabind_f_test_fixture <- function(methods = "HK", common = FALSE, n.by = 2L) {
  dat <- data.frame(
    event = rep(seq(10, 50, 10), n.by) + rep(20 * (seq_len(n.by) - 1), each = 5),
    total = 100, group = rep(LETTERS[seq_len(n.by)], each = 5))
  metaprop(event, total, data = dat, method = "GLMM", method.random.ci = methods,
             common = common, subgroup = group, tau.common = TRUE)
}

metabind_f_test_forest <- function(x, ...) {
  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off())
  meta::forest(x, layout = "meta", details = FALSE, ...)
}

subgroup_inference_fixture <- function(methods = "HK") {
  dat <- data.frame(event = c(10, 20, 30, 40, 50, 25, 35, 45, 55, 65),
                    total = 100, group = rep(c("A", "B"), each = 5))
  metaprop(event, total, data = dat, method = "GLMM",
           method.random.ci = methods, subgroup = group,
           common = TRUE, tau.common = TRUE)
}

subgroup_inference_print <- function(x, ...) {
  capture.output(print(x, details.methods = FALSE, digits.pval.Q = 4,
                       scientific.pval = FALSE, zero.pval = TRUE,
                        JAMA.pval = FALSE, ...))
}

subgroup_inference_forest <- function(x, ...) {
  grDevices::pdf(NULL, width = 16, height = 12)
  on.exit(grDevices::dev.off())
  grob <- grid::grid.grabExpr(
    forest(x, layout = "RevMan5", print.Q.subgroup = TRUE,
           label.test.subgroup.common = "INFERENCE_COMMON",
           label.test.subgroup.random = "INFERENCE_RANDOM",
           digits.pval.Q = 4, scientific.pval = FALSE,
            zero.pval = TRUE, JAMA.pval = FALSE, ...),
    width = 16, height = 12)
  labels <- function(g) {
    c(if (inherits(g, "text"))
        paste(deparse(g$label, width.cutoff = 500L), collapse = " "),
      unlist(lapply(g$children, labels), use.names = FALSE))
  }
  result <- labels(grob)
  result[grepl("INFERENCE_COMMON|INFERENCE_RANDOM", result)]
}

subgroup_difference_fixture <- function(methods = "HK", method = "GLMM") {
  dat <- data.frame(event = c(10, 20, 30, 40, 50, 30, 40, 50, 60, 70),
                    total = 100, group = rep(c("A", "B"), each = 5))
  metaprop(event, total, data = dat, method = method, method.random.ci = methods,
             common = FALSE, subgroup = group, tau.common = TRUE)
}

subgroup_difference_forest_label <- function(x, layout,
                                             tag = "SUBGROUP_DIFFERENCE", ...) {
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      args <- list(...)
      observed <<- lapply(args$cols$col.studlab$labels, `[[`, "label")
      do.call(draw, args)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 16, height = 12)
  on.exit(grDevices::dev.off())
  forest(x, layout = layout, label.test.subgroup.random = "SUBGROUP_DIFFERENCE",
         resid.hetlab = "RESIDUAL_HETEROGENEITY",
         print.Q.subgroup = TRUE, digits.Q = 2, digits.pval.Q = 4,
         scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE, ...)
  Filter(function(label) is.call(label) &&
           identical(label[[2]], tag), observed)
}

subgroup_residual_fixture <- function(methods = "classic", clustered = FALSE) {
  args <- list(
    TE = c(-0.4, -0.1, 0.2, 0.1, 0.6, 0.3, 0.5, 0.2, 0.9, 0.8, 1.1, 0.6),
    seTE = rep(0.2, 12), studlab = paste0("S", 1:12),
    subgroup = rep(c("A", "B"), each = 6), common = FALSE,
    tau.common = TRUE, method.tau = "REML", method.random.ci = methods)
  if (clustered)
    args$cluster <- rep(1:6, each = 2)
  do.call(metagen, args)
}

expect_subgroup_residual_print <- function(x) {
  output <- capture.output(print(
    x, details.methods = FALSE, digits.Q = 2, digits.pval.Q = 4,
    scientific.pval = FALSE, zero.pval = TRUE, JAMA.pval = FALSE))
  residual <- which(output == "Test of residual heterogeneity:")
  expect_length(residual, 1)
  expected <- paste(meta:::formatN(x$Q.resid, digits = 2),
                     meta:::formatN(x$df.Q.resid, digits = 0),
                     meta:::formatPT(x$pval.Q.resid, digits = 4,
                                    scientific = FALSE, zero = TRUE, JAMA = FALSE))
  expect_identical(strsplit(trimws(output[residual + 2]), "[[:space:]]+")[[1]],
                   strsplit(trimws(expected), "[[:space:]]+")[[1]])
  expect_false(any(grepl("^Within groups", output)))
  expect_equal(sum(grepl("^Test for subgroup differences", output)),
                length(x$method.random.ci))
  summary.output <- capture.output(print(summary(x), details.methods = FALSE))
  expect_equal(sum(summary.output == "Test of residual heterogeneity:"), 1)
  expect_false(any(grepl("^Within groups", summary.output)))
}

three_level_forest_fixture <- function(methods = "HK", rho = 0.3,
                                        three.level = TRUE) {
  dat <- data.frame(
    TE = c(0.1, 0.4, 0.2, 0.8, 0.3, 0.6, 1.2, 1.5, 0.9, 1.8, 1.1, 1.6),
    seTE = rep(c(0.1, 0.2, 0.15), 4), cluster = rep(1:6, each = 2),
    group = rep(c("B", "A"), each = 6))
  metagen(TE, seTE, data = dat, sm = "MD",
          cluster = if (three.level) cluster else NULL, rho = rho,
          subgroup = group, tau.common = TRUE, method.tau.ci = "",
          method.random.ci = methods)
}

three_level_forest_capture <- function(x, fields, ...) {
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      args <- list(...)
      observed <<- list(values = mget(fields, envir = parent.frame()), args = args)
      do.call(draw, args)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 20, height = 16)
  on.exit(grDevices::dev.off(), add = TRUE)
  forest(x, ...)
  observed
}

metabind_component_fields <- c(
  "tau2", "lower.tau2", "upper.tau2", "se.tau2",
  "tau", "lower.tau", "upper.tau")

metabind_component_fixture <- function(subgroup = FALSE) {
  dat <- data.frame(
    TE = c(-0.4, 0.1, 0.2, 0.7, -0.2, 0.4, 0.8, 1.1, 0.3, 0.9, 1.2, 0.6),
    seTE = rep(c(0.15, 0.2), 6), cluster = rep(1:6, each = 2),
    group = rep(c("A", "B"), each = 6))
  args <- list(TE = dat$TE, seTE = dat$seTE, data = dat, sm = "MD",
               common = FALSE, random = TRUE, tau.common = FALSE,
               method.random.ci = "HK", prediction = TRUE, method.predict = "V")
  if (subgroup)
    args$subgroup <- dat$group
  two <- do.call(metagen, args)
  args$cluster <- dat$cluster
  three <- do.call(metagen, args)
  sources <- list(two, three)
  # Name even scalar fields so extraction cannot silently discard names.
  for (i in seq_along(sources))
    for (field in metabind_component_fields) {
      value <- as.numeric(sources[[i]][[field]])
      sources[[i]][[field]] <- setNames(value, paste0("component", seq_along(value)))
    }
  sources
}

interval_slot_fixture <- function(labels = FALSE, shared = TRUE) {
  ma <- metagen(c(-0.4, 0.1, 0.7, 0.3, 0.9, 1.2, 0.2, 0.6),
                 rep(c(0.15, 0.2), 4), sm = "MD",
                 subgroup = rep(c("A", "B"), each = 4),
                 method.tau = "DL", method.random.ci = c("HK", "HK"),
                 adhoc.hakn.ci = "se", prediction = TRUE,
                 method.predict = c("V", "V"))
  if (labels) {
    ma$text.random <- c("First CI", "Second CI")
    ma$text.predict <- c("First PI", "Second PI")
  }
  x <- metabind(ma, ma, keepmeta = TRUE)
  if (!shared)
    x$samedata <- FALSE
  x
}

interval_slot_draw <- function(x, ...) {
  observed <- NULL
  draw <- meta:::forest.meta
  testthat::local_mocked_bindings(forest.meta = function(...) {
    observed <<- list(...)
    do.call(draw, observed)
  }, .package = "meta")
  grDevices::pdf(NULL, width = 24, height = 18)
  on.exit(grDevices::dev.off())
  result <- forest(x, layout = "meta", ...)
  list(args = observed, result = result)
}

metabind_subgroup_fixture <- function() {
  metabin(1:9, 101:109, 2:10, 201:209,
           subgroup = rep(c("A", "B", "C"), each = 3),
           method.random.ci = rep("HK", 3),
           adhoc.hakn.ci = c("", "se", "ci"),
           prediction = TRUE,
           method.predict = c("V", "HK", "HK", "HK-PR"),
           adhoc.hakn.pi = c("", "", "se", "se"))
}

metabind_repeated_interval_fixture <- function(labels = TRUE) {
  ma <- metabin(1:9, 101:109, 2:10, 201:209,
                  subgroup = rep(c("A", "B", "C"), each = 3),
                  method.random.ci = c("HK", "HK"), adhoc.hakn.ci = "se",
                  prediction = TRUE, method.predict = c("V", "V"),
                  adhoc.hakn.pi = "")
  if (labels) {
    ma$text.random <- c("First HK", "Second HK")
    ma$text.predict <- c("First PI", "Second PI")
  }
  metabind(ma, ma, keepmeta = TRUE)
}
