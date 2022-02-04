#' Smoking example
#' 
#' @description
#' Meta-analyses on the effect of smoking on mortality risk.
#' 
#' @details
#' Data have been reconstructed based on the famous Smoking and Health
#' Report to the Surgeon General (Bayne-Jones S et al., 1964). Data
#' sets can be used to evaluate the risk of smoking on overall
#' mortality (dataset \code{smoking}) and lung-cancer deaths (dataset
#' \code{lungcancer}), respectively.
#'
#' The person time is attributed such that the rate ratios are equal
#' to the reported mortality ratios implicitly assuming that the data
#' have arisen from a homogeneous age group; more detailed information
#' by age is not available from the report. Note, the group of
#' "non-smokers" actually consists of all participants except those
#' who are smokers of cigarettes only. Information on real non-smokers
#' is not available from the published Smoking and Health Report.
#' 
#' @name smoking
#' 
#' @aliases smoking lungcancer
#' 
#' @docType data
#' 
#' @format A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{study}}\tab study label \cr
#' \bold{\emph{participants}}\tab total number of participants \cr
#' \bold{\emph{d.smokers}}\tab number of deaths in smokers' group \cr
#' \bold{\emph{py.smokers}}\tab person years at risk in smokers' group
#'   \cr
#' \bold{\emph{d.nonsmokers}}\tab number of deaths in non-smokers'
#'   group \cr
#' \bold{\emph{py.nonsmokers}}\tab person years at risk in
#'   non-smokers' group
#' }
#' 
#' @seealso \code{\link{metainc}}
#' 
#' @source
#' Bayne-Jones S et al. (1964):
#' Smoking and Health: Report of the Advisory Committee to the Surgeon
#' General of the United States.
#' U-23 Department of Health, Education, and Welfare.
#' Public Health Service Publication No. 1103.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(smoking)
#' 
#' m1 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
#'   data = smoking, studlab = study)
#' print(m1, digits = 2)
#' 
#' data(lungcancer)
#' 
#' m2 <- metainc(d.smokers, py.smokers, d.nonsmokers, py.nonsmokers,
#'   data = lungcancer, studlab = study)
#' print(m2, digits = 2)


NULL
