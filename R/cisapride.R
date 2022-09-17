#' Cisapride in Non-Ulcer Dispepsia
#' 
#' @description
#' Meta-analysis on cisapride in non-ulcer dispepsia.
#' 
#' This meta-analysis is used as a data example in Hartung and Knapp
#' (2001).
#' 
#' @name cisapride
#' 
#' @docType data
#' 
#' @format A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study label \cr
#' \bold{\emph{event.cisa}}\tab number of events in cisapride group
#'   \cr
#' \bold{\emph{n.cisa}}\tab number of observations in cisapride group
#'   \cr
#' \bold{\emph{event.plac}}\tab number of events in placebo group \cr
#' \bold{\emph{n.plac}}\tab number of observations in placebo group
#' }
#' 
#' @seealso \code{\link{metabin}}
#' 
#' @source
#' Hartung J & Knapp G (2001):
#' A refined method for the meta-analysis of controlled clinical
#' trials with binary outcome.
#' \emph{Statistics in Medicine},
#' \bold{20}, 3875--89
#' 
#' @keywords datasets
#' 
#' @examples
#' data(cisapride)
#' 
#' m.or <- metabin(event.cisa, n.cisa, event.plac, n.plac,
#'   data = cisapride, sm = "OR",
#'   method = "Inverse", method.tau = "DL",
#'   studlab = study, method.incr = "all")
#' 
#' m.or.hk <- update(m.or, method.random.ci = "HK")
#' m.rr <- update(m.or, sm = "RR")
#' m.rr.hk <- update(m.or, sm = "RR", method.random.ci = "HK")
#'
#' vars.common <- c("TE.common", "lower.common", "upper.common")
#' vars.random <- c("TE.random", "lower.random", "upper.random")
#' #
#' res.common.or <- as.data.frame(m.or[vars.common])
#' names(res.common.or) <- vars.random
#' #
#' res.common.rr <- as.data.frame(m.rr[vars.common])
#' names(res.common.rr) <- vars.random
#' 
#' # Results for log risk ratio - see Table VII in Hartung and Knapp (2001) 
#' #
#' res.rr <- rbind(res.common.rr,
#'   as.data.frame(m.rr[vars.random]),
#'   as.data.frame(m.rr.hk[vars.random]))
#' #
#' row.names(res.rr) <- c("CE", "RE", "RE (HaKn)")
#' names(res.rr) <- c("Log risk ratio", "CI lower", "CI upper")
#' #
#' res.rr
#' 
#' 
#' # Results for log odds ratio (Table VII in Hartung and Knapp 2001) 
#' #
#' res.or <- rbind(res.common.or,
#'   as.data.frame(m.or[vars.random]),
#'   as.data.frame(m.or.hk[vars.random]))
#' #
#' row.names(res.or) <- c("CE", "RE", "RE (HaKn)")
#' names(res.or) <- c("Log odds ratio", "CI lower", "CI upper")
#' #
#' res.or


NULL
