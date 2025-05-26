###################################################
## Functions for: GPAbin: unifying visualizations of multiple mputations for 
## missing value
## Communications in Statistics: Simulation and Computation
## Authors: J Nienkemper-Swanepoel, NJ le Roux & S Gardner-Lubbe
## Corresponding author: J Nienkemper-Swanepoel
## 2023
###################################################


###################################################################################
biplFig <- function (CLPs, Zs, Lvls=NULL, Z.col="grey37", CLP.col="forestgreen", Z.pch=1, CLP.pch=17,title="") 
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function constructs a biplot after MCA
  ###################################################################################
  #Arguments
  #"CLPs" category level points (standard coordinates for variables)
  #"Zs" sample principal coordinates
  #"Lvls" names of the CLPs
  #"Z.col", "CLP.col" are the colour specifications for samples and CLPs
  #"Z.pch", "CLP.pch" are the plotting character specifications for samples and CLPs
  #"title" title of the figure
  ###################################################################################
  #Value
  #Returns a biplot.
  ###################################################################################
  ###################################################################################
  dev.new()
  par(pty="s")
  plot(rbind(CLPs[,1:2],Zs[,1:2]),pch="",xaxt="n",yaxt="n",xlab="",ylab="",main=title)
  points(Zs,pch=Z.pch,col=Z.col)
  points(CLPs,pch=CLP.pch,col=CLP.col)
  
  is.null(Lvls)
  {
    text(CLPs,cex=0.7,label=rownames(CLPs),pos=3)
  }
  !is.null(Lvls)
  {
    text(CLPs,cex=0.7,label=Lvls,pos=3)
  }
}

###################################################################################
###################################################################################

CLPred <- function (comp.lvls, datIN, CLPs, Zs, nsamples, pvar) 
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function predicts category levels for an MCA based biplot using the distances #between samples and CLPs
  ###################################################################################
  #Arguments
  #"comp.lvls" is the number of levels per variable in the original simulated data #set
  #"datIN" is input data set (could contain missing values)
  #"CLPs" and "Zs" contains the coordinate matrices of the GPAbin biplot
  #"nsamples" is the number of samples (rows) in the data sets
  #"pvar" is the number of variables (columns) in the data set
  ###################################################################################
  #Value
  #"Pred" a final predicted categorical data set
  ###################################################################################
  #Required functions
  #Auxiliary functions: is.integer0(), df2fact(), delCL(), FormatDimNam()
  ###################################################################################
  ###################################################################################
  I <- nsamples#number of samples
  
  comp.nam <-  levels(comp.lvls[,1])[comp.lvls[,1]]
  comp.nams <- unique(substr(comp.nam,1,1))
  NA.nams <- unique(substr(rownames(CLPs),1,1))
  finder <- which(is.na(match(comp.nams,NA.nams)))
  
  if (is.integer0(finder))
  {
    datIN <- datIN
  } else
  {
    datIN <- datIN[,-finder]
  }
  
  J <- ncol(datIN)
  X.pred <- as.data.frame(matrix(0,I,J))
  X.pred <- df2fact(X.pred)
  
  vec.lvl <- vector("numeric",J)
  for (j in 1:J)
  {
    vec.lvl[j] <- length(levels(delCL(datIN)[,j]))
  }
  cumlvl <- cumsum(vec.lvl)
  
  for (i in 1:I)
  {
    d <- as.matrix(as.matrix(dist(rbind(Zs[i,],CLPs)))[1,-1])#distance matrix between first row of Z.list[[i]][1,] and all rows of CLP.list[[i]] (distances for first sample over all variables)
    pluggin <- matrix(0,J,1)
    inds <- matrix(0,J,1)
    for(j in 1:J) 
    {
      if(j==1)
      {
        if(length(j:cumlvl[j])==1)
        {inds[j] <- 1
        } else
        {
          pluggin[j] <- min(d[j:cumlvl[j]])
          inds[j] <- which(d[j:cumlvl[j]]==min(d[j:cumlvl[j]]),arr.ind=TRUE)
        } 
      }else
      {
        if(length((cumlvl[j-1]+1):cumlvl[j])==1)
        {inds[j] <- 1
        }else
        {pluggin[j] <- min(d[(cumlvl[j-1]+1):cumlvl[j]])
        inds[j] <- which(d[(cumlvl[j-1]+1):cumlvl[j]]==min(d[(cumlvl[j-1]+1):cumlvl[j]]),arr.ind=TRUE)
        }
      } 
      levels(X.pred[,j]) <- levels(delCL(datIN)[,j])
      place <- inds[j]
      X.pred[i,j] <- levels(X.pred[,j])[place]
    }
  }
  inds.NNA <- which(!is.na(datIN), arr.ind=T)
  X.pred[inds.NNA] <- datIN[inds.NNA]#replaces non missing with original categories
  X.pred <- FormatDimNam(X.pred)
  
  return(Pred=X.pred)
}

###################################################################################
###################################################################################

compMeas <- function (Target=comp.dat , Testee=imp.dat, dim=c("All", "2D"),pred.dat1=pred.comp, pred.dat2=pred.imp,nsamples, pvar) 
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function calculates the following measures of comparison between two #configurations
  #Sim. pct, PS, CC, AMB, MB, RMSB
  ###################################################################################
  #Arguments
  #"Target" is the target configuration
  #"Testee" is the testee configuration
  #"dim" to compare the configurations in 2D or the maximum available ("All") dims
  #"pred.dat1" and "pred.dat2" are the predicted responses obtained from the CLPred function
  #"nsamples" number of samples in original data set
  #"pvar" number of variables original data set
  ###################################################################################
  #Value
  #Returns a data frame with the measures of comparison as stated above.
  ###################################################################################
  ###################################################################################
  counter <- 0
  Target <- as.matrix(Target)
  Testee <- as.matrix(Testee)
  
  nCLTar <- nrow(Target)
  nCLTes <- nrow(Testee)
  
  Tarnam <- rownames(Target)
  Tesnam <- rownames(Testee)
  
  #finding the CLs that occur in both Target and Testee and deleting the CLs that do #not appear in both in order to obtain a one-to-one comparison 
  rem <- which(is.na(match(Tarnam,Tesnam)))
  if(is.integer0(rem))
  {
    Target <- Target
    counter <- counter+1#counts the matched cases
  } else {Target<- Target[-rem,]}
  
  if(dim=="All")
  {
    pY <- ncol(Target)
    pX <- ncol(Testee)
    #the maximum number of common columns to use
    colUse <- min(pY,pX)
    Target <- Target[,1:colUse]
    Testee <- Testee[,1:colUse]
  } else
    if (dim=="2D")
    {
      Target <- Target[,1:2]
      Testee <- Testee[,1:2]
    }
  
  OPA <- myOPA (Y.in=Target,X.in=Testee,centring=TRUE)
  PS <- OPA[[2]]
  CC <- sum(dist(Testee) * dist(Target))/(sqrt(sum(dist(Testee)^2)) * sqrt(sum(dist(Target)^2)))
  RMSB <- ((sum(sum((Target-Testee)^2)))/length(Testee))^(0.5)
  MB <- (sum(sum((Target-Testee)^1)))/length(Testee)
  AMB <- (sum(sum(abs(Target-Testee))))/length(Testee)
  FitSS <- OPA[[5]]
  ResSS <- OPA[[3]]
  TotSS <- OPA[[4]]
  
  match.count <- sum(mapply(as.character,pred.dat1)==mapply(as.character,pred.dat2))
  
  sim.pct <- match.count/(pvar*nsamples)
  
  REStable <- data.frame(c(sim.pct, PS, CC, AMB, MB, RMSB))
  colnames(REStable)<- c("Measure of comparison")
  rownames(REStable)<- c("Sim.pct","PS", "CC", "AMB", "MB", "RMSB")
  print(REStable)
  return(REStable)
}

###################################################################################
###################################################################################

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

###################################################################################
###################################################################################

df2fact <- function (in.df, ordered=FALSE) 
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function transforms each column of a data.frame into unordered factor #variables (ordered=F) or ordered factor variables (ordered=T)
  ###################################################################################
  #Arugments
  #"in.df" a data.frame containing factors
  #"ordered" if TRUE (T) factor levels are ordered, if FALSE (F) factor levels are #unordered ###################################################################################
  #Value
  #Returns a transformed data.frame according to "ordered" argument.
  ###################################################################################
  ###################################################################################
  out.df <- factor(in.df[,1],ordered=ordered)
  
  for(i in 2:ncol(in.df)) out.df <- data.frame(out.df, factor(in.df[,i], ordered=ordered))
  
  colnames(out.df) <- colnames(in.df)
  rownames(out.df) <- rownames(in.df)
  return(out.df)
}

###################################################################################
###################################################################################

FormatDat <- function (datNA) 
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function is used to change the names of rows, columns and factor levels for #consistent display of results
  ###################################################################################
  #Arguments
  #"datNA" a multivariate categorical data set containing missing values
  ###################################################################################
  #Value
  #Returns a formatted "datNA"
  ###################################################################################
  ###################################################################################
  nn <- nrow(datNA)
  pp <- ncol(datNA)
  let.vec <- c(LETTERS[1:8],LETTERS[10:26])#excluding I
  lvl.vec <- c("one", "two", "thr","fou", "fiv", "six", "sev","eig", "nin", "ten","ele","twe","thi","frt","fivt")
  
  for (i in 1:nn) {rownames(datNA) <- paste('s',1:nn,sep='')}
  for (j in 1:pp) {colnames(datNA) <- let.vec[1:pp]
  numb <- length(levels(datNA[,j]))
  levels(datNA[,j]) <-lvl.vec[1:numb]
  }
  datNA
}

###################################################################################
###################################################################################

FormatDimNam <- function(dat) 
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function changes row- and column names, but not factor levels
  ###################################################################################
  #Arguments
  #"dat" a multivariate categorical data set
  ###################################################################################
  #Value
  #Returns formatted "dat"
  ###################################################################################
  ###################################################################################
  nn <- nrow(dat)
  pp <- ncol(dat)
  let.vec <- c(LETTERS[1:8],LETTERS[10:26])       #excluding I
  
  for (i in 1:nn) {rownames(dat) <- paste('s',1:nn,sep='')}
  for (j in 1:pp)
  {
    colnames(dat) <- let.vec[1:pp]
  }
  return(dat)
}

###################################################################################
###################################################################################

FormatImpList <- function (mylist)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function adds names to elements of the list containing the imputations
  ###################################################################################
  #Arguments
  #"mylist" is a list containing MIs
  ###################################################################################
  #Value
  #Returns formatted "mylist"
  ###################################################################################
  ###################################################################################
  imp <- length(mylist)
  names(mylist) <- paste('Imputation',1:imp,sep='.')
  return(mylist)
}

###################################################################################
###################################################################################

GPA <- function(Xk, G.target=NULL, iter=500, eps=0.001)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function contains the OPA function to compare two configurations and the GPA #function for multiple configuration comparisons
  ####################################################################################Arguments
  #"Xk" argument is a list containing the testee configurations which is updated on #each iteration
  #"G.target" argument is the target configuration. If not specified the centroid #configuration will be used as the target
  #"iter" is the number of iterations allowed before convergence
  #"eps" is the threshold value for convergence of the alogrithm
  ###################################################################################
  #Value
  #"Xk.F" is a list containing the updated testee configurations
  #"sk.F" is a vector containing the final scaling factors
  #"Qk.F" is a list containing the final rotation matrices
  #"Gmat" is the final target configuration
  #"sum.sq" is the final minimised sum of squared distance
  ###################################################################################
  ###################################################################################
  OPA <- function(X.mat, Z.mat)
  {
    svd.zx <- svd(t(Z.mat) %*% X.mat)
    svd.zx[[3]] %*% t(svd.zx[[2]])
  }
  K <- length(Xk)
  n <- nrow(Xk[[1]])
  p <- ncol(Xk[[1]])
  means <- t(sapply(Xk, function(X)apply(X, 2, mean)))
  Xk.scale <- sapply (Xk, scale, simplify=F)
  Xk.F <- sapply (1:K, function(k, Xk) as.matrix (Xk[[k]]), Xk=Xk.scale, simplify=F)
  Qk <- sapply(1:K, function(k, p) return(diag(p)), p=p, simplify = F)
  sk <- rep(1, K)
  sk.F <- sk
  Qk.F <- sapply(1:K, function(k, Qk) Qk[[k]], Qk = Qk, simplify = F)
  tel <- 0
  sum.sq.old <- Inf
  
  repeat
  {tel <- tel + 1
  if(tel > iter)
    stop(paste("Maximum number of specified iterations reached! Increase iter \n",II))
  Xk.F <- sapply(1:K, function(k, Xk, sk, Qk) sk[k] * Xk[[k]] %*% Qk[[k]], Xk = Xk.F, sk = sk, Qk = Qk, simplify = F)
  if (is.null(G.target))
  {Gmat <- Xk.F[[1]]
  for(k in 2:K)  Gmat <- Gmat + Xk.F[[k]]
  Gmat <- Gmat/K
  }
  else Gmat <- G.target 
  Qk <- sapply(1:K, function(k, Xind, Gmat) OPA(Xind[[k]], Gmat), Xind = Xk.F, Gmat = Gmat, simplify = F)
  Qk.F <- sapply(1:K, function(k, Qk, QF) QF[[k]] %*% Qk[[k]], Qk = Qk, QF = Qk.F, simplify = F)
  Smat <- matrix(0, ncol = K, nrow = K)
  for(i in 1:K)
    for(j in i:K) {
      Smat[i, j] <- sum(diag(t(Qk[[i]]) %*% t(Xk.F[[i]]) %*% Xk.F[[j]] %*% Qk[[j]]))
      if(i != j) Smat[j, i] <- Smat[i, j]
    }
  Smat.min.half <- diag(1/sqrt(diag(Smat)))
  swd <- svd(Smat.min.half %*% Smat %*% Smat.min.half)
  sk <- Smat.min.half %*% swd[[2]][, 1] * sqrt(K)
  sk.factor <- sum(sk)/K
  sk <- sk / sk.factor
  if(sk[1] < 0) sk <- -1 * sk
  
  sum.sq <- sum(sapply(1:K, function(k, Xk, Gmat)sum(diag((Xk[[k]] - Gmat) %*% t(Xk[[k]] - Gmat))), Xk = Xk.F, Gmat = Gmat))
  #cat("iter", tel, "sum.sq: ", sum.sq, "\n")
  if((sum.sq.old - sum.sq) < eps) break
  sum.sq.old <- sum.sq
  }
  list(Xk.F=Xk.F, sk.F=sk.F, Qk.F=Qk.F, Gmat=Gmat, sum.sq=sum.sq)
}

###################################################################################
###################################################################################

GPAbin <- function(CLP.list,Z.list,G.target=NULL)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function combines multiple configurations obtained from the output of the #MI.impute() function
  ####################################################################################Arguments
  #"CLP.list" argument is the list contains the CLPs of the multiple imputations
  #"Z.list" argument is the list contains the sample points of the multiple #imputations
  #"G.target" argument is the target configuration. If not specified the centroid #configuration will be used as the target
  ###################################################################################
  #Value
  #"Z.GPAbin" is the sample coordinates for the GPAbin biplot
  #"CLP.GPAbin" is the CLPs for the GPAbin biplot
  #"Z.GPA.list" is a list containing the sample coordinates for each MI after GPA
  #"CLP.GPA.list" is a list containing the CLPs for each MI after GPA
  ###################################################################################
  #Required functions
  #GPA()
  ###################################################################################
  ###################################################################################
  M <- length(CLP.list)
  GPA.out <- GPA(CLP.list,G.target=NULL)
  G.target <- GPA.out[[4]]
  Q.list <- GPA.out[[3]]
  s.list <- GPA.out[[2]]
  CLP.GPA.list <- GPA.out[[1]]
  
  Z.scal <- vector("list",M)
  for (m in 1:M)
  {
    Z.scal[[m]] <- Z.list[[m]]*s.list[[m]]
  }
  
  Z.GPA.list <- vector("list",M)
  for (m in 1:M)
  {
    Z.GPA.list[[m]] <- Z.scal[[m]]%*%Q.list[[m]]
  }
  
  Z.GPAbin <- Reduce("+",Z.GPA.list)/length(Z.GPA.list)
  CLP.GPAbin <- Reduce("+",CLP.GPA.list)/length(CLP.GPA.list)
  
  return(list(Z.GPAbin=Z.GPAbin, CLP.GPAbin=CLP.GPAbin, CLP.GPA.list=CLP.GPA.list, Z.GPA.list=Z.GPA.list,G.target=G.target))
}

###################################################################################
###################################################################################

indcol <- function(col.vec)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function constructs the dummy variables for the columns of an indicator #matrix per variable
  #This function is used in combination with the indmat() function
  ###################################################################################
  #Arguments
  #"col.vec" is a particular column from a multivariate categorical data set
  ###################################################################################
  #Value
  #Returns the dummy coded variables for a particular column from a multivariate #categorical data set
  ###################################################################################
  ###################################################################################
  elements <- levels(factor(col.vec))
  Y <- matrix(0, nrow = length(col.vec), ncol = length(elements))
  dimnames(Y) <- list(NULL, paste(elements))
  for(i in 1:length(elements))
  {
    Y[col.vec == elements[i], i] <- 1
  }
  return(Y)
}

indmat <- function(dat)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function constructs an indicator matrix
  ###################################################################################
  #Arguments
  #"dat" is a multivariate categorical data set
  ###################################################################################
  #Value 
  #Returns an indicator matrix
  ###################################################################################
  #Required functions
  #Auxiliary function: indcol()
  ###################################################################################
  ###################################################################################
  cols <- ncol(dat)
  samps <- nrow(dat)
  out <- matrix(1,samps,1)
  for (k in 1:cols)
  {
    ncol <- indcol(dat[,k])
    out <- cbind(out,ncol)
  }
  out <- out[,-1]
  return(out)
}

###################################################################################
###################################################################################

is.integer0 <- function(x)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function detects integers with zero length
  #This function is used as a precaution to eliminate possible errors
  ###################################################################################
  #Argument
  #"x" scalar value
  ###################################################################################
  #Value
  #Returns a TRUE or FALSE
  ###################################################################################
  ###################################################################################
  is.integer(x) && length(x) == 0L
}

###################################################################################
###################################################################################

MIimpute <- function(datNA=NULL, seed=123, imps=10)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function contains the function calls for MI using the MIMCA algorithm and #performing MCA on the MIs
  ###################################################################################
  #Arguments
  #"datNA" is a categorical data set with missing data entries
  #"seed" fixes the random seed in order to replicate results
  #"imps" specifies the number of multiple imputations
  ###################################################################################
  #Value
  #"CLP.list" is a list containing the CLPs of the MCA biplots of the MI data sets
  #"Z.list" is a list containing the sample coordinates of the MCA biplots of the MI #data sets
  #"datNA" is the input missing data
  #"Imp.list" is a list containing the MI categorical data sets
  ###################################################################################
  #Required packages
  #ca and FactoMineR
  ###################################################################################
  #Required functions
  #myestim_ncpMCA() and myimputeMCA()
  #Auxiliary functions: FormatDat(), indmat(), indcol(), FormatImpList(), rmOneCL()
  ###################################################################################
  ###################################################################################
  require(ca)
  require(FactoMineR)
  
  datNA <- as.data.frame(datNA)
  datNA <- FormatDat(datNA) #formatting row and column names
  set.seed(seed)
  colZscal <- ncol(indmat(datNA))#determining the number of columns in the #indicator matrix
  pvar <- ncol(datNA)#number of variables in missing data set
  ncp.max <- colZscal-pvar-1#maximum number of dimensions available for MCA #solution
  m <- imps
  #error handling of myestim_npcMCA()
  #stop() has been updated in myimputeMCA() to return the proposed number of #dimensions to retain and is used to run the function again without manual input
  
  nd <-try(myestim_ncpMCA(datNA,method="Regularized",method.cv="Kfold",ncp.min=0, ncp.max=ncp.max,verbose=FALSE, nbsim = 100,seed=seed)[[1]],silent=T)
  if(inherits(nd,"try-error")) nd<-as.numeric(conditionMessage(attr(nd,"condition")))
  
  set.seed(seed)
  
  if(is.na(nd))
  {
    #change seed if nd cannot be estimated
    nd <- try(myestim_ncpMCA(datNA,method="Regularized",method.cv="Kfold",ncp.min=0, ncp.max=ncp.max,verbose=FALSE, nbsim = 100,seed=(saad+12345))[[1]],silent=F)
    if(inherits(nd,"try-error")) nd<-as.numeric(conditionMessage(attr(nd,"condition")))
  }
  set.seed(seed)
  
  MI.output <- try(myMIMCA(datNA, ncp=nd, nboot=m, verbose=FALSE, seed=seed), silent=F)
  if(inherits(MI.output,"try-error")) MI.output <- myMIMCA(datNA, ncp=(as.numeric(conditionMessage(attr(MI.output,"condition")))),nboot=m, verbose=F, seed=seed)
  
  set.seed(seed)
  
  Imp.list <- MI.output[[1]]#list of imputed data sets
  Imp.list <- FormatImpList(Imp.list)#preparing colnames and rownames
  Imp.list <- rmOneCL(Imp.list)#preparation for MCA
  
  Z.list <- vector("list",m)
  CLP.list <- vector("list",m)
  
  for (imp in 1:m)
  {
    dat <- Imp.list[[imp]]
    out <- mjca(dat,lambda="indicator")
    nam <- out[[6]]
    Z.list[[imp]] <- out[[16]]
    CLPs <- out[[23]]
    rownames(CLPs) <- nam
    CLP.list[[imp]] <- CLPs
  }
  return(list(CLP.list=CLP.list,Z.list=Z.list,datNA=datNA, Imp.list=Imp.list,ncp=nd))
}

###################################################################################
###################################################################################

myestim_ncpMCA <- function (don, ncp.min = 0, ncp.max = 5, method = c("Regularized", "EM"), method.cv = c("Kfold", "loo"), nbsim = 100, pNA = 0.05, threshold = 1e-03, verbose = TRUE,seed=seed) 
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function includes minor changes to the original estim_ncpMCA() function #available in the missMDA R package
  #Code starting and ending with #change# indicates changes made to original functions
  ###################################################################################
  #Changes made to arguments
  #"threshold" changed to 1e-0.3
  #"seed" argument added
  ###################################################################################
  #Value
  #Returns the number of dimensions to be used
  ###################################################################################
  ###################################################################################
  #change# addition of fixed seed
  set.seed(seed)
  #change#
  tab.disjonctif.NA <- function(tab)
  {
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) 
    {
      moda <- tab[, i]
      nom <- names(tab)[i]
      n <- length(moda)
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind <- (1:n) + n * (unclass(moda) - 1)
      indNA <- which(is.na(ind))
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA, ] <- NA
      
      if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), "n", "N", "y", "Y"))) 
        dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), sep = "."))
      else dimnames(x) <- list(row.names(tab), levels(moda))
      return(x)
    }
    if (ncol(tab) == 1) 
      res <- modalite.disjonctif(1)
    else {
      res <- lapply(1:ncol(tab), modalite.disjonctif)
      res <- as.matrix(data.frame(res, check.names = FALSE))
    }
    return(res)
  }
  
  prodna <- function(x, noNA,seed=NULL) 
  {
    #change# addition of fixed seed
    set.seed(seed)
    #change#
    n <- nrow(x)
    p <- ncol(x)
    NAloc <- rep(FALSE, n * p)
    NAloc[sample(n * p, floor(n * p * noNA))] <- TRUE
    x[matrix(NAloc, nrow = n, ncol = p)] <- NA
    return(x)
  }
  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"), several.ok = T)[1]
  method.cv <- match.arg(method.cv, c("loo", "Kfold", "kfold", "LOO"), several.ok = T)[1]
  method <- tolower(method)
  method.cv <- tolower(method.cv)
  auxi = NULL
  don <- droplevels(don)
  
  for (j in 1:ncol(don)) if (is.numeric(don[, j])) 
    auxi = c(auxi, colnames(don)[j])
  
  if (!is.null(auxi)) 
    stop(paste("\nAll variables are not categorical, the following ones are numeric: ", auxi))
  vrai.tab = tab.disjonctif.NA(don)
  
  if (method.cv == "kfold")
  {
    res = matrix(NA, ncp.max - ncp.min + 1, nbsim)
    
    if (verbose) 
      pb <- txtProgressBar(min = 1/nbsim * 100, max = 100, style = 3)
    for (sim in 1:nbsim) 
    {
      continue <- TRUE
      while (continue) {
        donNA <- prodna(don, pNA,seed)
        continue <- (sum(unlist(sapply(as.data.frame(donNA), 
                                       nlevels))) != sum(unlist(sapply(don, nlevels))))
      }
      #change#
      #Omit for loop, factor levels are ordered incorrectly
      #for (i in 1:ncol(don)) donNA[, i] = as.factor(as.character(donNA[,i]))
      #change#
      for (nbaxes in ncp.min:ncp.max)
      {
        tab.disj.comp <- myimputeMCA(as.data.frame(donNA), ncp = nbaxes, method = method, threshold = threshold)[[1]]
        if (sum(is.na(donNA)) != sum(is.na(don))) 
        {
          res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE)/(sum(is.na(tab.disjonctif.NA(donNA))) - sum(is.na(tab.disjonctif.NA(don))))
        }  else
        {
          #change#
          #no alternative provided for cases where the sum of missing values are equal in don #and donNA res is a matrix of NAs, therefore else statement was added.
          
          donNA <- prodna(don, pNA,seed+1)
          continue <- (sum(unlist(sapply(as.data.frame(donNA), nlevels))) != sum(unlist(sapply(don, nlevels))))
          #change#
          
          #change#
          #Omit for loop, since factor levels are ordered incorrectly
          #for (i in 1:ncol(don)) donNA[, i] = as.factor(as.character(donNA[,i]))
          #change#
          for (nbaxes in ncp.min:ncp.max)
          {
            tab.disj.comp <- myimputeMCA(as.data.frame(donNA), ncp = nbaxes, method = method, threshold = threshold)[[1]]
            if (sum(is.na(donNA)) != sum(is.na(don))) 
            {
              res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE)/(sum(is.na(tab.disjonctif.NA(donNA))) - sum(is.na(tab.disjonctif.NA(don))))
            }
          }
        }
      }
      if (verbose) 
        setTxtProgressBar(pb, sim/nbsim * 100)
    }
    
    if (verbose) 
      close(pb)
    crit = apply(res, 1, mean, na.rm = TRUE)
    
    names(crit) <- c(ncp.min:ncp.max)
    result = list(ncp = as.integer(which.min(crit) + ncp.min - 1),
                  criterion = crit)
    return(result)
  }
  #change#
  #code omitted for “loo” method in appendix, since not of interest for this study
}

###################################################################################
###################################################################################

myimputeMCA <- function (don, ncp = 2, method = c("Regularized", "EM"), row.w = NULL, coeff.ridge = 1, threshold = 1e-03, seed = NULL, maxiter = 1000) 
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function includes minor changes to the original imputeMCA () function #available in the missMDA R package
  #Code starting and ending with #change# indicates changes made to original functions
  ###################################################################################
  #Value
  #"tab.disj" is the indicator matrix with fuzzy values for the missing values #according to the proportion of response category levels per variable
  #"completeObs" is the final categorical data set after SI
  ###################################################################################
  #Changes made to arguments
  #"threshold" changed to 1e-0.3
  ###################################################################################
  #Required packages
  #FactoMineR
  ###################################################################################
  ###################################################################################
  moy.p <- function(V, poids)
  {
    res <- sum(V * poids, na.rm = TRUE)/sum(poids[!is.na(V)])
  }
  
  find.category <- function(X, tabdisj)
  {
    nbdummy <- rep(1, ncol(X))
    is.quali <- which(!unlist(lapply(X, is.numeric)))
    
    #change#
    #To only consider the available observed category levels
    for (i in is.quali)
    {
      X[,i] <- droplevels(X[,i],exclude=NA)
    }
    #change#
    nbdummy[is.quali] <- unlist(lapply(X[, is.quali, drop = FALSE], nlevels))
    vec = c(0, cumsum(nbdummy))
    Xres <- X
    
    for (i in is.quali)
    {
      #change#
      #Original function did not foresee variables with only one category level observed. #This could occur when the percentage of missing values increase
      if (length((vec[i] + 1):vec[i + 1])==1)
      {
        temp <- as.factor(levels(X[,i]))
      } else
      {
        #change#
        temp <- as.factor(levels(X[, i])[apply(tabdisj[,(vec[i] + 1):vec[i + 1]], 1, which.max)])
      }
      Xres[, i] <- factor(temp, levels(X[, is.quali][,i]))
    }
    return(Xres)
  }
  
  tab.disjonctif.NA <- function(tab) 
  {
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) 
    {
      moda <- tab[, i]
      nom <- names(tab)[i]
      n <- length(moda)
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind <- (1:n) + n * (unclass(moda) - 1)
      indNA <- which(is.na(ind))
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA, ] <- NA
      if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), 
                                                     "n", "N", "y", "Y"))) 
        dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), sep = "."))
      else dimnames(x) <- list(row.names(tab), levels(moda))
      return(x)
    }
    if (ncol(tab) == 1) 
      res <- modalite.disjonctif(1)
    else  {
      res <- lapply(1:ncol(tab), modalite.disjonctif)
      res <- as.matrix(data.frame(res, check.names = FALSE))
    }
    return(res)
  }
  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"), 
                      several.ok = T)[1]
  method <- tolower(method)
  don <- droplevels(don)
  if (is.null(row.w)) 
    row.w <- rep(1/nrow(don), nrow(don))
  if (ncp == 0) 
    return(list(tab.disj = tab.disjonctif.prop(don, NULL, row.w = row.w), completeObs = find.category(don, tab.disjonctif.prop(don, NULL, row.w = row.w))))
  tab.disj.NA <- tab.disjonctif.NA(don)
  hidden <- which(is.na(tab.disj.NA))
  tab.disj.comp <- tab.disjonctif.prop(don, seed, row.w = row.w)
  tab.disj.rec.old <- tab.disj.comp
  continue <- TRUE
  nbiter <- 0
  while (continue) {
    nbiter <- nbiter + 1
    M <- apply(tab.disj.comp, 2, moy.p, row.w)/ncol(don)
    if (any(M < 0)) 
      #change#
      #stop statement updated to only provide the proposed number of dimensions
      #this enables the use of call function to run the function again with the proposed #number of dimensions without requiring manual feedback
      
      stop(ncp-1)
    
    #stop(paste("The algorithm fails to converge. Choose a number of components (ncp) #less or equal than " ,ncp - 1, " or a number of iterations (maxiter) less or equal #than ",maxiter - 1, sep = ""))
    #change#
    Z <- t(t(tab.disj.comp)/apply(tab.disj.comp, 2, moy.p, row.w))
    Z <- t(t(Z) - apply(Z, 2, moy.p, row.w))
    Zscale <- t(t(Z) * sqrt(M))
    svd.Zscale <- svd.triplet(Zscale, row.w = row.w, ncp = ncp)
    moyeig <- 0
    if (nrow(don) > (ncol(Zscale) - ncol(don))) 
      moyeig <- mean(svd.Zscale[[1]][-c(1:ncp, (ncol(Zscale)-ncol(don) + 1):ncol(Zscale))]^2)
    else moyeig <- mean(svd.Zscale[[1]][-c(1:ncp)]^2)
    moyeig <- min(moyeig * coeff.ridge, svd.Zscale[[1]][ncp + 1]^2)
    if (method == "em") 
      moyeig <- 0
    eig.shrunk <- ((svd.Zscale[[1]][1:ncp]^2 - moyeig)/svd.Zscale[[1]][1:ncp])
    if (ncp == 1) 
      rec <- tcrossprod(svd.Zscale[[2]][, 1] * eig.shrunk, svd.Zscale[[3]][, 1])
    else rec <- tcrossprod(t(t(svd.Zscale[[2]][, 1:ncp, drop = FALSE]) * 
                               eig.shrunk), svd.Zscale[[3]][, 1:ncp, drop = FALSE])
    tab.disj.rec <- t(t(rec)/sqrt(M)) + matrix(1, nrow(rec), ncol(rec))
    tab.disj.rec <- t(t(tab.disj.rec) * apply(tab.disj.comp, 2, moy.p, row.w))
    diff <- tab.disj.rec - tab.disj.rec.old
    diff[hidden] <- 0
    relch <- sum(diff^2 * row.w)
    tab.disj.rec.old <- tab.disj.rec
    tab.disj.comp[hidden] <- tab.disj.rec[hidden]
    continue = (relch > threshold) & (nbiter < maxiter)
  }
  
  tab <- find.category(don, tab.disj.comp)
  
  return(list(tab.disj = tab.disj.comp, completeObs = tab))
}

###################################################################################
###################################################################################

myMIMCA <- function (X, nboot = 10, ncp, coeff.ridge = 1, threshold = 1e-03, 
                     maxiter = 1000, verbose = FALSE, seed=NULL)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function includes minor changes to the original MIMCA () function available #in the missMDA R package
  #Code starting and ending with #change# indicates changes made to original functions
  ###################################################################################
  #Changes made to arguments
  #"threshold" changed to 1e-0.3
  ###################################################################################
  #Value
  #"X.imp" is a list containing the MI categorical data sets
  #"tab.disj" is the indicator matrix with fuzzy values for the missing values #according to the proportion of response category levels per variable
  #A list of of the input arguments are given:
  #"X", "nboot", "ncp", "coeff.ridge", "threshold", "seed", "maxiter", "tab.disj"
  ###################################################################################
  #Required packages
  #FactoMineR
  ###################################################################################
  #Required functions
  #myimputeMCA() function
  ###################################################################################
  ###################################################################################
  imputeMCA.print <- function(don, ncp, method = c("Regularized", "EM"), row.w = NULL, coeff.ridge = 1, threshold = 1e-03, seed = NULL, maxiter = 1000, verbose, printm)
  {
    
    if (verbose)
    {
      cat(paste(printm, "...", sep = ""))
    }
    res <- myimputeMCA(don = don, ncp = ncp, method = method, row.w = row.w, coeff.ridge = coeff.ridge, threshold = threshold, seed = seed, maxiter = maxiter)
    
    return(res)
  }
  normtdc <- function(tab.disj, data.na) {
    tdc <- tab.disj
    tdc[tdc < 0] <- 0
    tdc[tdc > 1] <- 1
    col.suppr <- cumsum(sapply(data.na, function(x) {nlevels(x)}))
    tdc <- t(apply(tdc, 1, FUN = function(x, col.suppr) {
      if (sum(x[1:col.suppr[1]]) != 1) {
        x[1:col.suppr[1]] <- x[1:col.suppr[1]]/sum(x[1:col.suppr[1]])
      }
      for (i in 2:length(col.suppr))
      {
        x[(col.suppr[i - 1] + 1):(col.suppr[i])] <- x[(col.suppr[i-1] + 1):(col.suppr[i])] /sum(x[(col.suppr[i-1] + 1):col.suppr[i]])
      }
      return(x)
    }, col.suppr = col.suppr))
    return(tdc)
  }
  
  draw <- function(tabdisj, Don, Don.na) {
    nbdummy <- rep(1, ncol(Don))
    is.quali <- which(!unlist(lapply(Don, is.numeric)))
    
    #change# edit to handle only one CL
    for (i in is.quali)
    {
      Don[,i] <- droplevels(Don[,i],exclude=NA)
    }
    #change#
    nbdummy[is.quali] <- unlist(lapply(Don[, is.quali, drop = FALSE],nlevels))
    vec = c(0, cumsum(nbdummy))
    Donres <- Don
    
    #change#
    for (i in is.quali) {
      if (length((vec[i] + 1):vec[i + 1])==1)
      {
        temp <- as.factor(levels(Don[,i]))
      } else
        #change#
        
      {
        Donres[, i] <- as.factor(levels(Don[, i])[apply(tabdisj[,(vec[i] + 1):vec[i+1]], 1, function(x) {
          sample(1:length(x), size = 1, prob = x)
        })])
      }
      Donres[, i] <- factor(Donres[, i], levels(Don[, is.quali][, i]))
    }
    return(don.imp = Donres)
  }
  
  temp <- if (coeff.ridge == 1) {
    "regularized"
  }
  else if (coeff.ridge == 0) {"EM"}
  else {
    paste("coeff.ridge=", coeff.ridge)
  }
  if (verbose)
  {
    cat("Multiple Imputation using", temp, "MCA using", nboot, "imputed arrays", "\n")
  }
  n <- nrow(X)
  Boot <- matrix(sample(1:n, size = nboot * n, replace = T), n, nboot)
  Weight <- matrix(1/(n * 1000), n, nboot, dimnames = list(1:n, paste("nboot=", 1:nboot, sep = "")))
  Boot.table <- apply(Boot, 2, table)
  for (i in 1:nboot) {
    Weight[names(Boot.table[[i]]), i] <- Boot.table[[i]]
  }
  
  Weight <- sweep(Weight, 2, STATS = colSums(Weight), FUN = "/")
  Weight <- as.data.frame(Weight)
  
  # res.imp <- mapply(Weight, FUN = imputeMCA.print, MoreArgs = list(don = X, 
  #     ncp = ncp, coeff.ridge = coeff.ridge, method = "Regularized", 
  #     threshold = threshold, maxiter = maxiter, verbose = verbose), 
  #     printm = as.character(1:nboot), SIMPLIFY = FALSE)
  
  #change# calling myimputeMCA()
  res.imp <- mapply(Weight, FUN = myimputeMCA, MoreArgs = list(don = X, ncp = ncp, coeff.ridge = coeff.ridge, method = "Regularized", threshold = threshold, maxiter = maxiter), SIMPLIFY = FALSE)
  #change#
  
  tdc.imp <- lapply(res.imp, "[[", "tab.disj")
  res.comp <- lapply(res.imp, "[[", "completeObs")
  tdc.norm <- mapply(FUN = normtdc, tab.disj = tdc.imp, data.na = res.comp, SIMPLIFY = F)
  
  X.imp <- mapply(FUN = draw, tabdisj = tdc.norm, Don = res.comp, MoreArgs = list(Don.na = X), SIMPLIFY = F)
  
  res <- list(res.MI = X.imp, res.imputeMCA = myimputeMCA(X, ncp = ncp, coeff.ridge = coeff.ridge, threshold = threshold, seed = NULL, maxiter = maxiter)[[1]], call = list(X = X, nboot = nboot, ncp = ncp, coeff.ridge = coeff.ridge, threshold = threshold, seed = NULL, maxiter = maxiter, tab.disj = array(unlist(tdc.imp), dim = c(nrow(tdc.imp[[1]]), ncol(tdc.imp[[1]]), length(tdc.imp)))))
  class(res) <- c("MIMCA", "list")
  if (verbose) {
    cat("\ndone!\n")
  }
  return(res)
}

###################################################################################
###################################################################################

myOPA <- function(Y.in = comp.CLP, X.in = NULL, centring = TRUE) 
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function performs Orthogonal Procrustes Analysis on centred data
  ###################################################################################
  #Arguments
  #"Y.in" target configuration (complete simulated configuration when available)
  #"X.in" testee configuration
  #"centring" by centring the data before OPA, translation step is redundant
  ###################################################################################
  #Value
  #"X.new" is the updated testee configuration
  #"ProStat" is the Procrustes Statistic
  #"Res.SS" is the residual sum of squares (Gower & Dijksterhuis 2004)
  #"Tot.SS" is the total sum of squares (Gower & Dijksterhuis 2004)
  #"Fitted.SS" is the fitted sum of squares (Gower & Dijksterhuis 2004)
  ###################################################################################
  ###################################################################################
  n.Y <- nrow(Y.in)
  p.Y <- ncol(Y.in)
  n.X <- nrow(X.in)
  p.X <- ncol(X.in)
  
  X.in <- as.matrix(X.in)
  Y.in <- as.matrix(Y.in)
  
  if(!centring)
  {
    X.in <- X.in
    Y.in <- Y.in
  }
  else
  {
    X.in <- scale(X.in,T,F)
    #centre=T, scale=F results are similar to Cox and Cox, Gower and #Dijkersthuis, Borg and Groenen
    Y.in <- scale(Y.in,T,F)
  }
  #transformations
  C.mat <- t(Y.in)%*%X.in
  svd.C <- svd(C.mat)
  A.mat <- svd.C[[3]]%*%t(svd.C[[2]])
  s.fact <- sum(diag(t(Y.in)%*%X.in%*%A.mat))/sum(diag(t(X.in)%*%X.in))
  #Gower and Dijksterhuis P32
  b.fact <- as.vector(1/n.Y * t(Y.in - s.fact * X.in %*% A.mat)%*%rep(1,n.Y))
  
  X.new <- b.fact + s.fact*X.in%*%A.mat
  
  Res.SS <- sum(diag(t(((s.fact*X.in%*%A.mat)-Y.in))%*%((s.fact*X.in%*%A.mat)-Y.in)))
  Tot.SS <- s.fact^2*sum(diag(t(X.in)%*%X.in))+sum(diag(t(Y.in)%*%Y.in))
  Fitted.SS <- 2*s.fact*sum(diag(svd.C[[1]]))
  PS <- Res.SS/sum(diag(t(Y.in)%*%Y.in))
  
  return(list(X.new=X.new, ProStat=PS, Res.SS=Res.SS, Tot.SS=Tot.SS, Fitted.SS=Fitted.SS))
}

###################################################################################
###################################################################################

rmOneCL <- function(inlist)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function removes columns with only one CL before applying MCA
  ###################################################################################
  #Arguments
  #"inlist" list of factor objects
  ###################################################################################
  #Value
  #Returns the input list, now with removed columns
  ###################################################################################
  ###################################################################################
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