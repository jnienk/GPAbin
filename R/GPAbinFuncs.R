#' Deleting empty category levels
#' 
#' This function removes unobserved category levels before applying multivariate techniques.
#'
#' @param data is an object of the factor class
#'
#' @return
#' \item{data}{Returns the input object, now with empty category levels removed}
#' @export
#'
delCL <- function(data)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function removes empty factor levels
  ###################################################################################
  #Arguments
  #"data" is an object of the factor class
  ###################################################################################
  #Value
  #"Returns the input object, now with empty category levels removed
  ###################################################################################
  ###################################################################################
  for (j in 1:ncol(data))
  {
    col <- data[,j]
    nl <- levels(col)
    
    for (i in 1:length(nl))
    {
      col <- droplevels(col,exclude=NA)
    }
    data[,j] <- col
  }
  return(data)
}

