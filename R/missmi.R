### JNS
### GPAbin package functions
### ========================================================

.onAttach <- function(libname,pkgname){
  packageStartupMessage("Welcome to GPAbin! \nThis package is used to construct biplots for multiple imputed data sets. \nRun ?missmi for more information")
}

#' First step before constructing unified biplots
#' 
#' @description
#' This function produces a list of elements to be used when producing a GPAbin biplot.
#'
#' @param data input data frame or list
#'
#' @returns
#' \item{X}{The processed data}
#' \item{imputations}{Number of multiple imputations applied}
#' \item{n}{The number of samples}
#' \item{p}{The number of variables}
#' \item{miss_pct}{Percentage of missing values}
#'  
#' @export
#' 
#' @examples
#' data(missdat)
#' missbp <- missmi(missdat)
#' data(implist)
#' missbp <- missmi(implist)
#' 
missmi <- function(data)
{
  #think about separating continuous and categorical data for applications of PCA GPAbin
    
    if(!is.data.frame(data))
    {
      X <- lapply(data, apply, 2, as.factor) #ensuring all columns are factors
      m <- length(X)
      n <- nrow(X[[1]]) #assuming all lists will be the same size
      p <- ncol(X[[1]]) #assuming all lists will be the same size
      miss_pct <- NULL
    }
    else
    {
      X <- data
      n <- nrow(X)
      p <- ncol(X)
      miss_pct <- colSums(is.na(X))/n*100
      miss_over <- sum(miss_pct)/p
      m <- NULL
      if((miss_over==0) && is.null(m)) stop("The data is complete.")
    }

  missbp <- list(X = X, m = m, n = n, p=p, miss_pct = miss_pct)
  class(missbp) <- "missmi"
  
  missbp
}

###################################################################################
###################################################################################
#' Generic print function for objects of class missmi
#' @description
#' This function is used to print output when the missmi biplot object is created.
#' 
#' @param x an object of class \code{missmi}.
#' @param ... additional arguments.
#'
#' @return This function will not produce a return value, it is called for side effects.
#'
#' @export
#' @examples
#' data(missdat)
#' missbp <- missmi(missdat)
#' data(implist)
#' missbp <- missmi(implist)
#' print(missbp)
#' 
print.missmi <- function (x, ...)
{
  if(!is.data.frame(x$X))
  {
    print(paste("There are", x$m, "imputations / variations of your data available."))
  } else
  {
    print(paste("Missing data with a total of", round(mean(x$miss_pct),0), "% missing values."))
  }
}
