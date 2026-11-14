forest_named_regression_plot <- function(x, ...) {
  observed <- NULL
  draw <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest_meta_internal = function(...) {
      observed <<- list(...)
      do.call(draw, observed)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 18, height = 14)
  on.exit(grDevices::dev.off())
  result <- forest(x, ...)
  list(result = result, args = observed)
}

test_that("Cairo devices work as names and functions", {
  skip_if_not(capabilities("cairo"))
  ma <- metagen(1:3, 3:1)
  files <- c(tempfile(fileext = ".pdf"), tempfile(fileext = ".pdf"),
              tempfile(fileext = ".ps"), tempfile(fileext = ".ps"))
  on.exit(unlink(files))
  devices <- list("cairo_pdf", grDevices::cairo_pdf,
                   "cairo_ps", grDevices::cairo_ps)
  headers <- rep(c("%PDF-", "%!PS-"), each = 2)
  results <- vector("list", length(devices))
  #
  for (i in seq_along(devices)) {
    results[[i]] <- forest(ma, filename = files[i], device = devices[[i]])
    expect_true(file.exists(files[i]))
    expect_gt(file.info(files[i])$size, 0)
    expect_identical(rawToChar(readBin(files[i], "raw", n = 5)), headers[i])
  }
  for (i in c(1, 3)) {
    expect_equal(results[[i]]$width, results[[i + 1]]$width)
    expect_equal(results[[i]]$height, results[[i + 1]]$height)
  }
  expect_error(meta:::validate_device("unknown_device", files[1]),
                "Unknown graphics device")
})

test_that("character hetstat preserves complete model-row expressions", {
  datas <- data.frame(a_n = c(60, 70, 100), a_mean = c(8, 7, 10),
                      a_sd = c(11, 11, 13), study = c("s1", "s2", "s3"))
  ma <- metamean(a_n, a_mean, a_sd, study, data = datas)
  grouped <- metamean(rep(datas$a_n, 2), rep(datas$a_mean, 2),
                       rep(datas$a_sd, 2), subgroup = rep(c("A", "B"), each = 3))
  for (x in list(ma, grouped)) {
    for (placement in c("common", "fixed", "random")) {
      drawn <- forest_named_regression_plot(x, hetstat = placement)
      row <- if (placement == "random") 3 else 2
      label <- drawn$args$cols$col.studlab$labels[[row]]$label
      expect_type(label, "expression")
      expect_length(label, 1)
      expect_identical(label[[1]][[1]], quote(paste))
      expect_match(paste(deparse(label), collapse = " "), "Heterogeneity")
      expect_match(paste(deparse(label), collapse = " "), "italic\\(I\\)\\^2")
    }
  }
})

test_that("deprecated study hetstat warns and falls back to logical TRUE", {
  ma <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), method.tau = "DL")
  expected <- forest_named_regression_plot(ma, hetstat = TRUE)
  expect_warning(drawn <- forest_named_regression_plot(ma, hetstat = "study"),
                 paste0("Use of argument 'hetstat = \"study\"' is deprecated; ",
                        "using 'hetstat = TRUE' instead."),
                 fixed = TRUE)
  expect_identical(lapply(drawn$args$cols$col.studlab$labels, `[[`, "label"),
                   lapply(expected$args$cols$col.studlab$labels, `[[`, "label"))
  expect_identical(drawn$args$cols$col.studlab$rows,
                   expected$args$cols$col.studlab$rows)
})

test_that("omitted studies are removed from additional columns before sorting", {
  dat <- data.frame(te = c(NA, 0.2, 0.3), se = rep(0.1, 3),
                     extra = c("DROP", "KEEP2", "KEEP3"))
  ma <- metagen(te, se, data = dat, studlab = paste0("S", 1:3), method.tau = "DL")
  # A single source column must remain a data frame through selection and sorting.
  ma$data <- dat["extra"]
  drawn <- forest_named_regression_plot(ma, allstudies = FALSE, sortvar = -TE,
                                        leftcols = c("studlab", "extra", "lower", "effect", "ci"))
  expect_identical(drawn$result$studlab, c("S3", "S2"))
  labels <- lapply(drawn$args$cols$col.extra$labels, `[[`, "label")
  labels <- unname(unlist(Filter(function(z) is.character(z) &&
    length(z) == 1 && grepl("^(DROP|KEEP)", z), labels)))
  expect_identical(labels, c("KEEP3", "KEEP2"))
  expect_length(drawn$args$cols$col.extra$labels, length(drawn$args$yS))
  labels <- lapply(drawn$args$cols$col.lower$labels, `[[`, "label")
  labels <- unname(unlist(Filter(function(z) is.character(z) &&
    length(z) == 1 && grepl("^-?[0-9]", trimws(z)), labels)))
  expect_equal(trimws(labels), meta:::formatN(ma$lower[c(3, 2)], digits = 2))
})

test_that("single-column RoB data follow study removal and sorting", {
  ma <- metagen(c(NA, 0.2, 0.3), rep(0.1, 3), studlab = paste0("S", 1:3),
                 method.tau = "DL")
  risk <- rob(c("Low risk of bias", "High risk of bias", "Unclear risk of bias"),
                studlab = paste0("S", 1:3), tool = "rob1", domains = "Domain")
  drawn <- forest_named_regression_plot(ma, rob = risk, allstudies = FALSE,
                                        sortvar = -TE, rightcols = "RoB.A")
  expect_s3_class(drawn$args$rob, "data.frame")
  expect_equal(dim(drawn$args$rob), c(2, 1))
  expect_equal(as.character(drawn$args$rob[[1]]),
               c("Unclear risk of bias", "High risk of bias"))
  expect_identical(drawn$result$studlab, c("S3", "S2"))
})

test_that("subgroup layout respects inherited and model-specific types", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3), method.tau = "DL")
  for (legacy in c(FALSE, TRUE)) {
    if (legacy)
      drawn <- forest_named_regression_plot(ma, layout = "subgroup", type.subgroup = "circle")
    else
      drawn <- forest_named_regression_plot(ma, layout = "subgroup", type = list(subgroup = "circle"))
    expect_equal(drawn$args$col.forest$type[4:7], rep("circle", 4))
  }
  drawn <- forest_named_regression_plot(ma, layout = "subgroup",
                                        type = list(subgroup = "circle", subgroup.common = "diamond"))
  expect_equal(drawn$args$col.forest$type[4:7], c("diamond", "diamond", "circle", "circle"))
  drawn <- forest_named_regression_plot(ma, layout = "subgroup")
  expect_equal(drawn$args$col.forest$type[4:7], rep("square", 4))
})

test_that("legacy subgroup row spacing matches named-list settings", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3), method.tau = "DL")
  modern <- forest_named_regression_plot(ma, addrow = list(subgroup = FALSE))
  legacy <- forest_named_regression_plot(ma, addrow.subgroups = FALSE)
  expect_identical(legacy$args$cols$col.studlab$rows, modern$args$cols$col.studlab$rows)
  spaced <- forest_named_regression_plot(ma, addrow = list(subgroup = TRUE))
  expect_gt(spaced$args$nrow, modern$args$nrow)
  expect_error(forest_named_regression_plot(ma, addrow = list(subgroup = FALSE),
                                            addrow.subgroups = TRUE), "cannot be combined")
  expect_error(forest_named_regression_plot(ma, addrow.subgroup = FALSE,
                                            addrow.subgroups = TRUE), "cannot be combined")
})

test_that("automatic limits include finite vector CID thresholds", {
  ma <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), method.tau = "DL")
  drawn <- forest_named_regression_plot(ma, cid = c(1, 2), xlim = NULL)
  expect_lte(drawn$result$xlim[1], -2)
  expect_gte(drawn$result$xlim[2], 2)
  single <- metamean(rep(20, 3), c(0.1, 0.3, 0.5), rep(0.2, 3),
                      method.tau = "DL", null.effect = 0)
  drawn <- forest_named_regression_plot(single, cid = c(1, 2))
  expect_lte(drawn$result$xlim[1], -2)
  expect_gte(drawn$result$xlim[2], 2)
})

test_that("ratio CID thresholds use the back-transformed null effect", {
  data(Olkin1995, package = "meta")
  ma <- metabin(ev.exp, n.exp, ev.cont, n.cont,
                 data = Olkin1995, subset = c(41, 47, 51, 59),
                 sm = "RR", method = "I", studlab = paste(author, year))
  drawn <- forest_named_regression_plot(
    ma, layout = "RevMan5", common = FALSE,
    cid = 0.75, fill = list(fill = "lightgray", cid = "white"))
  expect_equal(drawn$args$ref, 0)
  expect_equal(drawn$args$cid.below.null, log(0.75))
  expect_equal(drawn$args$cid.above.null, log(1 / 0.75))
  #
  thresholds <- c(0.25, 0.5, 0.75)
  drawn <- forest_named_regression_plot(
    ma, cid = list(below.null = thresholds,
                    above.null = 1 / rev(thresholds)))
  expect_equal(drawn$args$cid.below.null, log(thresholds))
  expect_equal(drawn$args$cid.above.null, log(1 / rev(thresholds)))
})

test_that("CID reference conversion preserves nonzero and explicit references", {
  ma <- metagen(log(c(1, 1.5, 2)), rep(0.1, 3), sm = "RR",
                  null.effect = log(1.2), method.tau = "DL")
  drawn <- forest_named_regression_plot(
    ma, cid = list(below.null = 1.1, above.null = 1.3))
  expect_equal(drawn$args$ref, log(1.2))
  #
  drawn <- forest_named_regression_plot(
    ma, ref = 1.4, cid = list(below.null = 1.3, above.null = 1.5))
  expect_equal(drawn$args$ref, log(1.4))
  #
  drawn <- forest_named_regression_plot(ma, backtransf = FALSE)
  expect_equal(drawn$args$ref, log(1.2))
})

test_that("single-outcome CID references respect stored null-effect scales", {
  ma <- metaprop(c(2, 4, 8), c(100, 200, 300), method = "Inverse",
                  sm = "PLOGIT", null.effect = 0.2, method.tau = "DL")
  generic <- metagen(ma$TE, ma$seTE, sm = "PLOGIT",
                       null.effect = qlogis(0.2), method.tau = "DL")
  for (x in list(ma, generic)) {
    drawn <- forest_named_regression_plot(
      x, cid = list(below.null = 0.1, above.null = 0.3))
    expect_equal(drawn$args$ref, 0.2)
  }
})

test_that("explicit CID fill takes precedence over equivalence fill", {
  ma <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), method.tau = "DL")
  expect_warning(drawn <- forest_named_regression_plot(
    ma, cid = 0.2, fill = list(cid.below.null = "red", lower.equi = "blue")), "ignored")
  expect_identical(drawn$args$fill.cid.below.null[1], "red")
  expect_false("blue" %in% drawn$args$fill.cid.below.null)
})

test_that("stored component predictions and heterogeneity use original study counts", {
  for (singleton in c(FALSE, TRUE)) {
    if (singleton) {
      TE <- c(-0.4, 0.8)
      seTE <- c(0.1, 0.3)
      subgroup <- c("A", "B")
    }
    else {
      TE <- c(-0.4, 0.8, 0.5, -0.3)
      seTE <- c(0.1, 0.3, 0.2, 0.2)
      subgroup <- c("A", "A", "A", "B")
    }
    ma <- metagen(TE, seTE, subgroup = subgroup, common = FALSE, random = FALSE,
                   prediction = TRUE, method.predict = "V", method.tau = "DL")
    combined <- metabind(ma, update(ma, method.tau = "REML"),
                          name = c("DL component", "REML component"))
    for (allstudies in c(FALSE, TRUE)) {
      drawn <- forest_named_regression_plot(combined, allstudies = allstudies)
      labels <- lapply(drawn$args$cols$col.studlab$labels, `[[`, "label")
      heterogeneity <- vapply(labels, function(z)
        grepl("Heterogeneity", paste(deparse(z), collapse = " "), fixed = TRUE), FALSE)
      expect_equal(sum(heterogeneity), 2)
      expect_true(all(!is.na(drawn$args$cols$col.studlab$rows[heterogeneity])))
      expect_gte(sum(is.finite(drawn$args$col.forest$low)), 2)
      if (!allstudies)
        expect_identical(drawn$result$studlab,
                         if (singleton) character() else c("A", "A"))
    }
  }
})

test_that("risk-only settings retain the data and caption in the settings list", {
  ma <- metagen(c(0.1, 0.3, 0.5), rep(0.1, 3), studlab = paste0("S", 1:3),
                 method.tau = "DL")
  risk <- rob(c("Low risk of bias", "High risk of bias", "Unclear risk of bias"),
                studlab = paste0("S", 1:3), tool = "rob1", domains = "Domain")
  drawn <- forest_named_regression_plot(ma, rob = list(rob = risk, only = TRUE,
                                                       legend = FALSE, text = "Bias"))
  expect_true(drawn$args$RoB.available)
  labels <- lapply(drawn$args$cols$col.studlab$labels, `[[`, "label")
  expect_false(any(vapply(labels, function(z)
    grepl("Risk of bias legend", paste(deparse(z), collapse = " "), fixed = TRUE), FALSE)))
  expect_identical(drawn$result$rightcols, "RoB.A")
  expect_identical(drawn$args$col.rob$labels[[1]]$label, "Bias")
  drawn <- forest_named_regression_plot(ma, rob = FALSE)
  expect_false(drawn$args$RoB.available)
  expect_null(drawn$args$rob)
})

test_that("named-list settings render with legacy file height calculation", {
  ma <- metagen(seq(0.1, 0.6, 0.1), rep(0.1, 6),
                 subgroup = rep(c("A", "B"), each = 3), method.tau = "DL")
  file <- tempfile(fileext = ".pdf")
  withr::defer(unlink(file))
  expect_silent(forest(ma, filename = file, rows.gr = 12, width = 12,
                        label = list(left = "Favours treatment", right = "Favours comparator"),
                        fs = list(main = 14), ff = list(common = "italic"),
                        col = list(square = "navy"), addrow = list(subgroup = FALSE),
                        cid = list(cid = c(0.1, 0.2), pooled.only = TRUE),
                        test = list(overall = TRUE), details = TRUE))
  expect_true(file.exists(file))
  expect_gt(file.info(file)$size, 0)
})

test_that("dropping an empty subgroup retains the remaining heterogeneity counts", {
  ma <- metagen(c(NA, NA, 0.2, 0.3, 0.4), rep(0.1, 5),
                 subgroup = c("A", "A", "B", "B", "B"), method.tau = "DL")
  drawn <- forest_named_regression_plot(ma, allstudies = FALSE)
  labels <- vapply(drawn$args$cols$col.studlab$labels, function(grob)
    paste(deparse(grob$label), collapse = " "), "")
  visible <- !is.na(drawn$args$cols$col.studlab$rows)
  heterogeneity <- labels[visible & grepl("Heterogeneity", labels, fixed = TRUE)]
  expect_length(heterogeneity, 2)
  expect_false(any(grepl("not applicable", heterogeneity, fixed = TRUE)))
})
