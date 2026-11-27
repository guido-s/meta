test_that("metabin() applies single-study MH.exact fallback with double-zero studies", {
  for (sm in c("OR", "DOR", "RR", "VE")) {
    m1 <- metabin(
      event.e = 0, n.e = 40, event.c = 4, n.c = 40,
      sm = sm, method = "MH"
    )
    m2 <- metabin(
      event.e = c(0, 0, 0), n.e = c(40, 40, 40),
      event.c = c(0, 4, 0), n.c = c(40, 40, 40),
      sm = sm, method = "MH"
    )
    expect_false(m2$MH.exact)
    expect_equal(m2$k, 1)
    expect_equal(m2$TE.common, m1$TE.common)
    expect_equal(m2$seTE.common, m1$seTE.common)
    expect_equal(m2$lower.common, m1$lower.common)
    expect_equal(m2$upper.common, m1$upper.common)
    expect_equal(m2$pval.common, m1$pval.common)
  }

  # Double all-event studies are also excluded for OR when allstudies = FALSE
  m1_or <- metabin(
    event.e = 0, n.e = 40, event.c = 4, n.c = 40,
    sm = "OR", method = "MH"
  )
  m_ae <- metabin(
    event.e = c(40, 0, 40), n.e = c(40, 40, 40),
    event.c = c(40, 4, 40), n.c = c(40, 40, 40),
    sm = "OR", method = "MH"
  )
  expect_false(m_ae$MH.exact)
  expect_equal(m_ae$k, 1)
  expect_equal(m_ae$TE.common, m1_or$TE.common)
  expect_equal(m_ae$seTE.common, m1_or$seTE.common)

  # When double-zero studies are retained (allstudies = TRUE or sm = "RD"),
  # k > 1 so MH.exact stays TRUE
  m_all <- metabin(
    event.e = c(0, 0, 0), n.e = c(40, 40, 40),
    event.c = c(0, 4, 0), n.c = c(40, 40, 40),
    sm = "OR", method = "MH", allstudies = TRUE
  )
  expect_true(m_all$MH.exact)
  expect_equal(m_all$k, 3)

  m_rd <- metabin(
    event.e = c(0, 0, 0), n.e = c(40, 40, 40),
    event.c = c(0, 4, 0), n.c = c(40, 40, 40),
    sm = "RD", method = "MH"
  )
  expect_true(m_rd$MH.exact)
  expect_equal(m_rd$k, 3)
})
