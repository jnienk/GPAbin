#' Converting data frames to factors
#' 
#' Transforms each column of a data frame into factor variables.
#'
#' @param in.df a data.frame containing factors
#' @param ordered logical argument. If TRUE (T) factor levels are ordered, if FALSE (F) factor levels are unordered.
#'
#' @return  
#' \item{out.df}{Returns a transformed data frame according to "ordered" argument.}
#' 
#' @noRd
#'
df2fact <- function (in.df, ordered=FALSE)
{
  out.df <- factor(in.df[,1],ordered=ordered)
  
  for(i in 2:ncol(in.df)) out.df <- data.frame(out.df, factor(in.df[,i], ordered=ordered))
  
  colnames(out.df) <- colnames(in.df)
  rownames(out.df) <- rownames(in.df)
  return(out.df)
}
#########################################################################

#' Formatting Imputation Lists
#' 
#' Adding consistent names to list elements from imputations
#'
#' @param mylist list containing MIs
#'
#' @return
#' \item{mylist}{formatted `mylist`}
#' 
#' @noRd
#'
FormatImpList <- function (mylist)
{
  imp <- length(mylist)
  names(mylist) <- paste('Imputation',1:imp,sep='.')
  return(mylist)
}
###################################################################################

#' Removing empty category levels
#'
#' Removes columns with only one CL before applying dimension reduction
#' @param inlist list of factor objects
#'
#' @return
#' \item{inlist}{Returns the input list, now with removed columns}
#' 
#' @noRd
#'
rmOneCL <- function(inlist)
{
  M <- length(inlist)
  
  for(m in 1:M)
  {
    dat <- inlist[[m]]
    pvar <- ncol(dat)
    vecel <- vector("numeric",pvar)
    for (i in 1:pvar)
    {
      if(length(levels(dat[,i]))==1)
        vecel[i] <- i
    }
    if(sum(vecel)==0)
    {
      inlist[[m]] <- dat
    } else
    {
      vecel <- vecel[vecel != 0]
      inlist[[m]] <- dat[,-vecel]
    }
  }
  return(inlist)
}
###################################################################################
#' Removing empty category levels specifically for `jomo` imputation
#'
#' Removes columns with only one CL before applying MCA specifically for output from jomo
#' @param inlist list of factor objects
#'
#' @return
#' \item{inlist}{Returns the input list, now with removed columns}
#' 
#' @noRd
#'
rmOneCLjom <- function(inlist)
{
  M <- length(inlist)
  
  for(m in 1:M)
  {
    dat <- inlist[[m]][,1:5]
    pvar <- ncol(dat)
    vecel <- vector("numeric",pvar)
    for (i in 1:pvar)
    {
      if(length(levels(dat[,i]))==1)
        vecel[i] <- i
    }
    if(sum(vecel)==0)
    {
      inlist[[m]] <- dat
    } else
    {
      vecel <- vecel[vecel != 0]
      inlist[[m]] <- dat[,-vecel]
    }
  }
  return(inlist)
}
###################################################################################

#' Detects integers with zero length
#'
#' Used as a precaution to eliminate possible errors
#' @param x scalar value
#'
#' @return
#' TRUE or FALSE
#'
#' @noRd
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

###################################################################################

#' Removing empty factor levels
#'
#' Removes empty factor levels
#' @param data object of the factor class
#'
#' @return
#' \item{data} Returns the input object, now with empty category levels removed
#'
#' @noRd
#'
delCL <- function(data)
{
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