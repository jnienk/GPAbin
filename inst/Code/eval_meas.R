#fit measures
#continue work on this
#only to use if you have the complete data, as in a simulation study

evalMeas <- function(missbp, compdat, X.in = NULL, centring = TRUE)
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
  X.in <- rbind(missbp$GPAbin.Z, missbp$GPAbin.CLP)
  #X.in <- as.matrix(X.in) #testee GPAbin
  
  #DRT on compplete data
    temp <- ca::mjca(compdat,lambda="indicator")
    compZ <-temp[[16]]
    compCLP <- temp[[23]]
    complevels <- temp[[6]]

  Y.in <- rbind(compZ, compCLP)
  #Y.in <- as.matrix(Y.in) #target complete

  n.Y <- nrow(Y.in)
  p.Y <- ncol(Y.in)
  n.X <- nrow(X.in)
  p.X <- ncol(X.in)
  
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
  #Tot.SS <- s.fact^2*sum(diag(t(X.in)%*%X.in))+sum(diag(t(Y.in)%*%Y.in))
  #Fitted.SS <- 2*s.fact*sum(diag(svd.C[[1]]))
  PS <- Res.SS/sum(diag(t(Y.in)%*%Y.in))
  
  missbp$PS <- PS
  missbp$compZ <- compDRT$Z
  missbp$compCLP <- compDRT$CLP
  missbp$complevels <- compDRT$levels
  missbp$compdat <- compdat
  
  #return(list(X.new=X.new, ProStat=PS, Res.SS=Res.SS, Tot.SS=Tot.SS, Fitted.SS=Fitted.SS))
}
###################################################################################
CLPred <- function (missbp, CLPs=CLPs, Zs=Zs)
{
  ###################################################################################
  ###################################################################################
  #Information
  #This function predicts category levels for an MCA based biplot using the distances #between samples and CLPs
  ###################################################################################
  #Arguments
  #"mca.comp" MCA solution of original simulated case if available
  #"datIN" is input data set (could contain missing values)
  #"CLPs" and "Zs" contains the coordinate matrices of the GPAbin biplot
  #"nsamples" is the number of samples (rows) in the data sets
  #"pvar" is the number of variables (columns) in the data set
  ###################################################################################
  #Value
  #"Pred" a final predicted categorical data set
  ###################################################################################
  #Required functions
  #Auxiliary functions: is.integer0(), df2fact(), delCL()
  ###################################################################################
  ###################################################################################

  p <- missbp$p
  n <- missbp$n
  
  if(is.null(CLPs)) {
    CLPs <- missbp$CLP.GPAbin
  }else{
  }
  if(is.null(Zs)) {
    Zs <- missbp$Z.GPAbin
  }
  
  X.pred <- as.data.frame(matrix(0,n,p))
  X.pred <- df2fact(X.pred)
  
  vec.lvl <- vector("numeric",p)
  for (j in 1:p)
  {
    vec.lvl[j] <- length(levels(delCL(datIN)[,j]))
  }
  cumlvl <- cumsum(vec.lvl)
  
  for (i in 1:I)
  {
    d <- as.matrix(as.matrix(dist(rbind(Zs[i,],CLPs)))[1,-1]) #distance matrix between first row of Z.list[[i]][1,] and all rows of CLP.list[[i]] (distances for first sample over all variables)
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
  #X.pred <- FormatDimNam(X.pred)
  
  missbp$predCL <- missbp$X.pred
  
  missbp
}
