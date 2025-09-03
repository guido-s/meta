# for use with dplyr::group_map
pwgen_study <- function(input, studlab) {
  # Get rid of warning 'Undefined global functions or variables'
  treat1 <- NULL
  #
  studlab <- studlab$studlab

  # separate the reference arm from the other rows
  ref <- input[input$treat1 == input$treat2, ]
  stopifnot(nrow(ref) == 1)
  refSE <- ref$seTE
  ref <- ref$treat2
  rel <- input %>% filter(treat1 != ref) %>%
    mutate(studlab = studlab)
  stopifnot(nrow(rel) == nrow(input) - 1)
  stopifnot(nrow(rel) == length(unique(rel$treat1)))

  # calculate other treatment effects
  if (nrow(rel) > 1) {
    # generate the missing treatment pairs
    other <- combn(rel$treat1, 2)
    # for each pair, calculate TE and seTE
    other <- mapply(function(treat1, treat2) {
      r1 <- rel[rel$treat1 == treat1, ]
      r2 <- rel[rel$treat1 == treat2, ]
      data.frame(
        studlab,
        treat1,
        treat2,
        TE = r1$TE - r2$TE,
        seTE = sqrt(r1$seTE^2 + r2$seTE^2 - 2 * refSE^2))
    },
    other[1, ], other[2, ], SIMPLIFY = FALSE) %>% bind_rows()
    #
    bind_rows(rel, other)
  }
  else {
    rel
  }
}

# Expects one row per arm of each study.
#  - treat1 is the current study arm 
#  - treat2 is the study reference arm (e.g. control/placebo)
#  - if treat1 != treat2, TE is the treatment effect estimate, otherwise NA
#  - if treat1 != treat2, seTE is the standard error of the treatment effect
#    estimate, otherwise it is the standard error of the reference arm
#    (or the square root of the covariance between treatment effects) - optional
#    for two-arm studies.
pwgen <- function(studlab, treat1, treat2, TE, seTE) {
  data.frame(studlab, treat1, treat2, TE, seTE) %>%
    group_by(studlab) %>%
    group_map(pwgen_study) %>%
    bind_rows() %>%
    relocate(studlab) %>%
    as.data.frame
}

