metabind_empty_sm_fixture <- function(subgroup = FALSE, shared = TRUE) {
  dat <- data.frame(TE = c(0.1, 0.3, 0.5, 0.2, 0.4, 0.6),
                    seTE = rep(0.1, 6), group = rep(c("A", "B"), each = 3))
  args <- list(TE = dat$TE, seTE = dat$seTE, data = dat,
               common = TRUE, random = TRUE, prediction = TRUE,
               tau.common = FALSE, method.tau = "DL", method.predict = "V")
  if (subgroup)
    args$subgroup <- dat$group
  ma <- do.call(metagen, args)
  other <- if (shared) ma else update(ma, level.ma = 0.90)
  metabind(ma, other, name = c("First", "Second"), keepmeta = TRUE)
}

test_that("default empty summary measures print cleanly in every effect table", {
  for (scenario in c("plain", "shared", "nonshared")) {
    x <- metabind_empty_sm_fixture(scenario != "plain", scenario != "nonshared")
    before <- x
    if (scenario != "plain")
      expect_identical(meta:::prepare_metabind_subgroups(x)$shared.overall,
                       scenario == "shared")
    expect_identical(x$overall$sm, c("", ""))
    expect_silent(output <- capture.output(print(x, details.methods = FALSE)))
    expect_silent(summary.output <- capture.output(print(
      summary(x), details.methods = FALSE)))
    for (lines in list(output, summary.output)) {
      expect_false(any(grepl("\\.no_sm_provided|\\.\\.\\.[0-9]+", lines)))
      expect_true(all(c("Common effect model", "Random effects model", "Prediction interval") %in%
                        trimws(lines)))
      headers <- lines[grepl("CI|PI", lines) & !grepl("Prediction interval", lines)]
      expect_true(length(headers) > 0L)
      expect_false(any(grepl("\\bTE\\b|Effect", headers)))
      if (scenario == "shared")
        expect_true("Overall results" %in% lines)
    }
    if (scenario == "nonshared")
      expect_false("Overall results" %in% output)
    expect_identical(x, before)
    expect_true(all(vapply(x$meta, function(ma) identical(ma$sm, ""), logical(1))))
  }
})

test_that("empty effect headings are emitted as blanks, not dropped or repaired", {
  x <- metabind_empty_sm_fixture()
  flags <- list(k = FALSE, stat = FALSE, pval = FALSE, Q = FALSE, I2 = FALSE,
                tau2 = FALSE, tau = FALSE, H = FALSE)
  for (model in c("common", "random", "prediction")) {
    selection <- as.list(setNames(c(FALSE, FALSE, FALSE),
                                  c("common", "random", "prediction")))
    selection[[model]] <- TRUE
    expect_silent(output <- capture.output(do.call(
      print, c(list(x = x, print = flags, details.methods = FALSE), selection))))
    expect_false(any(grepl("\\.no_sm_provided|\\.\\.\\.[0-9]+|\\bTE\\b", output)))
    header <- output[grepl("95% (CI|PI)", output)]
    expect_length(header, 1L)
    expect_match(header, "^\\s+95% (CI|PI)")
    if (model != "prediction") {
      row <- output[grepl("^First\\s", output)]
      expect_length(row, 1L)
      effect <- regexpr("0.35", row, fixed = TRUE)[1]
      expect_gt(effect, 0L)
      expect_identical(substr(header, effect, effect + 3L), "    ")
    }
  }
  expect_identical(x$overall$sm, c("", ""))
})
