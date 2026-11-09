metabind_reference <- function(x, ...) {
  observed <- numeric()
  draw <- meta:::forest_meta_internal
  file <- tempfile(fileext = ".pdf")
  on.exit(unlink(file))

  testthat::with_mocked_bindings(
    forest(x, backtransf = FALSE, filename = file, ...),
    forest_meta_internal = function(...) {
      args <- list(...)
      observed <<- c(observed, args$ref)
      do.call(draw, args)
    },
    .package = "meta"
  )
  unique(observed)
}

test_that("single-outcome references are transformed once by forest.meta", {
  ma <- metaprop(c(2, 4, 8), c(100, 200, 300), method = "Inverse",
                 sm = "PLOGIT", null.effect = 0.2)
  mb <- metabind(ma, ma)
  expect_equal(metabind_reference(mb), qlogis(0.2))
  expect_equal(metabind_reference(mb, ref = 0.3), qlogis(0.3))

  generic <- metagen(ma$TE, ma$seTE, sm = "PLOGIT",
                      null.effect = qlogis(0.2))
  expect_equal(metabind_reference(metabind(generic, generic)), qlogis(0.2))
})

test_that("PFT references use unique harmonic means or the PAS fallback", {
  ma <- metaprop(c(2, 4, 8), c(100, 200, 300), method = "Inverse",
                 sm = "PFT", null.effect = 0.2)
  other <- metaprop(c(4, 8, 16), c(200, 400, 600), method = "Inverse",
                    sm = "PFT", null.effect = 0.2)
  expect_equal(metabind_reference(metabind(ma, ma)),
               meta:::transf(0.2, "PFT", n = ma$n.harmonic.mean))
  expect_equal(metabind_reference(metabind(ma, other)),
               meta:::transf(0.2, "PAS"))
})

test_that("IRFT references use unique harmonic times or a square-root fallback", {
  ma <- metarate(c(2, 4, 8), c(100, 200, 300), sm = "IRFT",
                 null.effect = 0.2)
  other <- metarate(c(4, 8, 16), c(200, 400, 600), sm = "IRFT",
                    null.effect = 0.2)
  expect_equal(metabind_reference(metabind(ma, ma)),
               meta:::transf(0.2, "IRFT", time = ma$t.harmonic.mean))
  expect_equal(metabind_reference(metabind(ma, other)), sqrt(0.2))
})
