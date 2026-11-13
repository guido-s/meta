forest_metabind_options <- function(x, ...) {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file, width = 18, height = 12)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  })
  meta::forest(x, layout = "meta", ...)
}

metagen_forest_options <- function(...) {
  metagen(log(c(1.2, 2, 3)), c(0.15, 0.2, 0.25),
          sm = "RR", method.tau = "DL", ...)
}

metabind_render_fixture <- function(method.tau = "DL", ...) {
  dat <- data.frame(
    event.e = c(6, 30, 12, 55, 8, 45, 18, 60),
    n.e = c(100, 110, 90, 120, 100, 110, 90, 120),
    event.c = c(20, 12, 25, 18, 24, 15, 30, 20),
    n.c = c(100, 110, 90, 120, 100, 110, 90, 120),
    first = rep(c("A", "B"), each = 4),
    second = rep(c("C", "D"), 4))
  ma <- metabin(event.e, n.e, event.c, n.c, data = dat,
                subgroup = first, sm = "RR", method = "Inverse",
                method.tau = "DL", common = TRUE, random = TRUE,
                prediction = TRUE, method.predict = "V", ...)
  other <- update(ma, subgroup = second, method.tau = method.tau)
  metabind(ma, other, name = c("First grouping", "Second grouping"),
           keepdata = TRUE, keepmeta = TRUE)
}

metabind_render_pdf <- function(x, ...) {
  observed <- NULL
  labels <- NULL
  columns <- NULL
  draw <- meta:::forest.meta
  draw.internal <- meta:::forest_meta_internal
  testthat::local_mocked_bindings(
    forest.meta = function(...) {
      observed <<- list(...)
      do.call(draw, observed)
    },
    forest_meta_internal = function(...) {
      args <- list(...)
      labels <<- lapply(args$cols$col.studlab$labels, `[[`, "label")
      columns <<- args$cols
      do.call(draw.internal, args)
    }, .package = "meta")
  grDevices::pdf(NULL, width = 24, height = 16)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- meta::forest(x, layout = "meta", ...)
  list(result = result, args = observed, labels = labels, columns = columns)
}

forest_transformation_columns <- function(x, ...) {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file, width = 18, height = 12)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  })
  meta::forest(x, ...)
}

metagen_transformation_columns <- function(...) {
  metagen(1:3, 3:1, sm = "MD", prediction = TRUE, ...)
}
