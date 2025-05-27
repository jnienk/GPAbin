#' Complete data example
#'
#' Simulated data example
#' 
#' \describe{
#'   \item{V1}{Variable 1}
#'   \item{V1}{Variable 2}
#'   \item{V1}{Variable 3}
#'   \item{V1}{Variable 4}
#'   \item{V1}{Variable 5}
#' }
#' 
#' @name compdat
#' @docType data
#' @format A data frame with 1000 rows and 5 columns.
#' @keywords datasets
#' @source Simulated data from a uniform distribution that is categorised into levels.
NULL

#' Missing data example
#'
#' `compdat` containing approximately 35% simulated missing values according to a missing at random (MAR) missing data mechanism
#' 
#' \describe{
#'   \item{V1}{Variable 1}
#'   \item{V1}{Variable 2}
#'   \item{V1}{Variable 3}
#'   \item{V1}{Variable 4}
#'   \item{V1}{Variable 5}
#' }
#' 
#' @name missdat
#' @docType data
#' @format A data frame with 1000 rows and 5 columns.
#' @keywords datasets
#' @source Simulated data from a uniform distribution that is categorised into levels.
NULL

#' List of multiple imputed data sets
#'
#' Five multiple imputations of `missdat`
#' 
#' \describe{
#'   \item{V1}{Variable 1}
#'   \item{V1}{Variable 2}
#'   \item{V1}{Variable 3}
#'   \item{V1}{Variable 4}
#'   \item{V1}{Variable 5}
#' }
#' 
#' @name implist
#' @docType data
#' @format List containing five multiple imputations of `missdat`. Each list item a data frame with 1000 rows and 5 columns.
#' @keywords datasets
#' @source simulated example data imputed with mice::mice(missdat, m=5, method="polyreg", maxit=10, remove.collinear=FALSE, printFlag = FALSE)
NULL

#' @rdname implist