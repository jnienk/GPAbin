### JNS
### GPAbin functions
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
missmi <- function(data)
{
  #think about separating continuous and categorical data for applications of PCA GPAbin
    
    if(is.data.frame(data[[1]]))
    {
      X <- lapply(data, apply, 2, as.factor) #ensuring all columns are factors
      
      m <- length(X)
      n <- nrow(X[[1]])
      p <- ncol(X[[1]])
      miss_pct <- NULL
      print(paste("There are", m, "imputations / variations of your data available."))
    }
    else
    {
      
      #X <- apply(data, 2, as.factor) #ensuring all columns are factors
      X <- data
      n <- nrow(X)
      p <- ncol(X)
      miss_pct <- colSums(is.na(X))/nrow(X)*100
      miss_over <- sum(miss_pct)/p
      m <- NULL
      if(sum(miss_pct==0) && is.null(m)) stop("The data is complete.")
      print(paste("Missing data with a total of ", round(miss_over,0), "% missing values."))
    }

  object <- list(X = X, imputations = m, n = n, p=p, miss_pct = miss_pct)
  class(object) <- "missmi"
  
  object
}