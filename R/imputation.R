#' Multiple imputation
#' 
#' Choose between four available multiple imputation strategies in `R`.
#'
#' @param missbp An object of class \code{missmi} obtained from preceding function \code{missmi()}.
#' @param imp.method Select one of four imputation methods: `MIMCA`, `jomo`, `miceB`, `mice`
#' @param m Number of multiple imputations
#'
#' @returns
#' The `missbp` object is appended with the following object:
#' \item{dataimp}{List of imputed data}
#' 
#' See also \code{\link[missMDA]{MIMCA}}, \code{\link[jomo]{jomo1cat}}, \code{\link[mi]{mi}} and \code{\link[mice]{mice}}.
#' 
#' @export
#'
#' @examples
#' \donttest{
#' data(missdat)
#' missbp <- missmi(missdat) |> impute(imp.method="miceB", m=5)}
#' 
impute <- function(missbp, imp.method=c("MIMCA","jomo","miceB","mice"), m=5)
{
  if(is.null(missbp$miss_pct)) stop("Don't apply an imputation method. Your data is already complete, continue to DRT().")
  
  n <- missbp$n
  p <- missbp$p
  missbp$m <- m
  data <- missbp$X
  
  if (imp.method=="MIMCA")
  {
    
    #Imputation
    ncp.out <- missMDA::estim_ncpMCA(data, method="Regularized", method.cv="Kfold") #regularisation step parameter
    MIMCA.imp <- missMDA::MIMCA(data, nboot=m, ncp=ncp.out[[1]]) #multiple imputation
    
    MIMCA.list <- MIMCA.imp[[1]] #list of m imputed data sets
    MIMCA.list <- FormatImpList(MIMCA.list) #preparing colnames and rownames
    MIMCA.list <- rmOneCL(MIMCA.list) #preparation for MCA
    
    missbp$dataimp <- MIMCA.list

  } else
    if (imp.method=="jomo")
    {
      ##Imputation
      cat.vec <- vector("numeric",ncol(data))
      for (j in 1:ncol(data))
      {
        cat.vec[j] <-length(levels(data[,j]))
      }
      
      jomo.out <- jomo::jomo1cat(data,cat.vec,nimp=m, output=0)
      jomo.list <- mitools::imputationList(split(jomo.out,jomo.out$Imputation)[-1])[[1]]
      
      samp.nams <- paste("s",1:n,sep="")
      jomo.list <- FormatImpList(jomo.list) #preparing colnames and rownames
      jomo.list <- rmOneCLjom(jomo.list) #preparation for MCA
      
      for (imp in 1:m)
      {
        rownames(jomo.list[[imp]]) <- samp.nams #uniform names across imputations
      }
      
      missbp$dataimp <- jomo.list
      
    } else
      if (imp.method=="miceB")
      {
        
        ##Imputation
        miceB.prep <- mi::missing_data.frame(data, subclass = "allcategorical")
        miceB.out <- mi::mi(miceB.prep, n.chains = m)
        miceB.out <- mi::complete(miceB.out, m = m)
        miceB.list <- lapply(miceB.out,'[',1:p)
        
        samp.nams <- paste("s",1:n,sep="")
        miceB.list <- rmOneCL(miceB.list) #preparation for MCA
        
        for (imp in 1:m)
        {
          rownames(miceB.list[[imp]]) <- samp.nams #uniform names across imputations
        }
        
        missbp$dataimp <- miceB.list
        
      } else
        if (imp.method=="mice")
        {
        
          ##Imputation
          mice.imp <- mice::mice(data,m=m,method="polyreg", maxit=10,remove.collinear=FALSE, printFlag = FALSE)
          mice.out <- mice::complete(mice.imp, "long")
          mice.list <- mitools::imputationList(split(mice.out,mice.out$.imp))[[1]]
          
          samp.nams <- paste("s",1:n,sep="")
          for (imp in 1:m)
          {
            mice.list[[imp]] <- mice.list[[imp]][,1:5] #preparing colnames and rownames
            rownames(mice.list[[imp]]) <- samp.nams #uniform names across imputations
          }
          mice.list <- rmOneCL(mice.list) #preparation for MCA

          missbp$dataimp <- mice.list
          
        } else {return("no imputation method provided")}
  
  missbp$m <- m
  
  missbp
}
###################################################################################