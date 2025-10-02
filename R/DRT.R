#' Dimension reduction
#' 
#' Multiple correspondence analysis is performed on the multiple imputed datasets
#'
#' @param missbp An object of class \code{missbp} obtained from preceding function \code{missmi()}
#' @param method Select a dimension reduction technique. In the current version `MCA` is available.
#'
#' @return
#' The `missbp` object is appended with the following objects:
#' \item{Z}{List of sample coordinates}
#' \item{CLP}{List of category level point coordinates}
#' \item{lvls}{List of category level names}
#' \item{m}{Number of multiple imputations}
#' 
#' See also \code{\link{missmi}} and \code{\link{impute}}.
#' @export
#'
#' @examples
#' data(implist)
#' missbp <- missmi(implist) |> DRT()
#' 
DRT <- function(missbp, method=c("MCA"))
{
  #in this version only multiple correspondence analysis
  
  #if completed data were provided
  if(is.null(missbp$miss_pct)) 
    {
    #continue with imputed data
    data <- missbp$X
  } else{
    data <- missbp$dataimp
    }
  
  m <- missbp$m
  
  Z.list <- vector("list", m)
  CLP.list <- vector("list", m)
  lvls.list <- vector("list", m)
  
  #ensure data is a factor
  for (i in 1:m)
  {
  data[[i]] <- df2fact(data[[i]])
  }
  #remove empty category levels prior to dimension reduction
  data <- rmOneCL(data)
  
  for (i in 1:m)
  {
  temp <- ca::mjca(data[[i]],lambda="indicator")
  Z.list[[i]] <-temp[[16]]
  CLP.list[[i]] <- temp[[23]]
  lvls.list[[i]] <- temp[[6]]
  }
  
  missbp$Z <- Z.list
  missbp$CLP <- CLP.list
  missbp$lvls <- lvls.list
  
  missbp
}
