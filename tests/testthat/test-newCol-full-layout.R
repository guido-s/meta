newCol_layout <- function(pooled, n.com = 1, n.ran = 2, n.prd = 2,
                          n.com.w = 2, n.ran.w = 4, n.prd.w = 4,
                          n.stat.w = 2, sel.prd = c(4, 5, 9, 10, 14, 15),
                          studies = c("study 1", "study 2")) {
  meta:::newCol("value", "Value", rob = NULL,
               data1 = data.frame(value = studies),
               data2 = NULL, datap = pooled,
               n.com = n.com, n.ran = n.ran, n.prd = n.prd,
               notavail = FALSE, lab.NA = "NA", big.mark = "",
               digits = 2, digits.tau2 = 2, digits.tau = 2, digits.I2 = 1,
               sel.prd = sel.prd,
               n.com.w = n.com.w, n.ran.w = n.ran.w, n.prd.w = n.prd.w,
               n.stat.w = n.stat.w)$format_var
}

test_that("full pooled layouts retain group-major model and prediction slots", {
  pooled <- c("overall common", "overall random 1", "overall random 2",
              "overall prediction 1", "overall prediction 2",
              "A common", "A random 1", "A random 2",
              "A prediction 1", "A prediction 2",
              "B common", "B random 1", "B random 2",
              "B prediction 1", "B prediction 2")
  expect_identical(newCol_layout(data.frame(value = pooled)),
                   c(pooled, "", "", "study 1", "study 2"))
  expect_identical(newCol_layout(data.frame(value = pooled), n.stat.w = 0),
                   c(pooled, "study 1", "study 2"))
})

test_that("full pooled layouts also support overall-only predictions", {
  expect_identical(
    newCol_layout(data.frame(value = c("common", "random", "prediction")),
                  n.ran = 1, n.prd = 1, n.com.w = 0, n.ran.w = 0,
                  n.prd.w = 0, n.stat.w = 0, sel.prd = 3),
    c("common", "random", "prediction", "study 1", "study 2"))
  expect_identical(
    newCol_layout(data.frame(value = "prediction"),
                  n.com = 0, n.ran = 0, n.prd = 1,
                  n.com.w = 0, n.ran.w = 0, n.prd.w = 0,
                  n.stat.w = 0, sel.prd = 1),
    c("prediction", "study 1", "study 2"))
})

test_that("full pooled layouts format numeric predictions and missing slots", {
  pooled <- seq_len(15) / 10
  pooled[c(9, 11)] <- NA_real_
  result <- newCol_layout(data.frame(value = pooled), studies = c(2, 3))
  expect_identical(result[c(4, 5, 10, 14, 15)],
                   c("0.40", "0.50", "1.00", "1.40", "1.50"))
  expect_identical(result[c(9, 11, 16, 17)], rep("", 4))
  expect_identical(result[18:19], c("2.00", "3.00"))
})

test_that("legacy pooled layouts retain their expansion and blank predictions", {
  expect_identical(newCol_layout(data.frame(value = "pooled")),
                   c(rep("pooled", 3), rep("", 14), "study 1", "study 2"))
  expect_identical(newCol_layout(data.frame(value = letters[1:3])),
                   c(letters[1:3], rep("", 14), "study 1", "study 2"))
  expected <- c(letters[1:3], "", "", letters[4:9], rep("", 6),
                "study 1", "study 2")
  expected[c(4, 5, 9, 10, 14, 15)] <- ""
  expect_identical(newCol_layout(data.frame(value = letters[1:9])), expected)
  expect_identical(newCol_layout(NULL),
                   c(rep("", 17), "study 1", "study 2"))
})
