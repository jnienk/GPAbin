###################################################################################
CLPpred <- function (CLPs=CLPs, Zs=Zs, p=p, n=n, lvls=lvls, datIN = datIN)
{
  #should be general for any given Z and CLP coordinates to predict CLs
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
  
  X.pred <- as.data.frame(matrix(0,n,p))
  X.pred <- df2fact(X.pred)
  
  vec.lvl <- vector("numeric",p)
  for (j in 1:p)
  {
    vec.lvl[j] <- length(levels(delCL(datIN)[,j]))
  }
  cumlvl <- cumsum(vec.lvl)
  
  for (i in 1:n)
  {
    d <- as.matrix(as.matrix(dist(rbind(Zs[i,],CLPs)))[1,-1]) #distance matrix between first row of Z.list[[i]][1,] and all rows of CLP.list[[i]] (distances for first sample over all variables)
    pluggin <- matrix(0,p,1)
    inds <- matrix(0,p,1)
    for(j in 1:p)
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
  
  return(predCL = X.pred)
  
}

#continue work on this
#only to use if you have the complete data, as in a simulation study


evalMeas <- function (missbp, compdat=NULL, dim=c("All", "2D"))
                      #Z.Target=Z.comp , CLP.Target=CLP.comp, Z.Testee=Z.imp, CLP.Testee=CLP.imp,
                   #   dim=c("All", "2D"),pred.dat1=pred.comp, pred.dat2=pred.imp,nsamples=n, pvar=vars)
{
  if(is.null(compdat)) stop("No complete data set available for comparison. \nThis function is only applicable for simulated data.")
  
  ###################################################################################
  ###################################################################################
  #Information
  #This function calculates the following measures of comparison between two #configurations
  #PS, SP, RPR, AMB, RMSB
  ###################################################################################
  #Arguments
  #"Z.Target" is the sample coordinates of target configuration
  #"CLP.Target" is the category level point coordinates of target configuration
  #"Z.Testee" is the sample coordinates of testee configuration
  #"CLP.Testee" is the category level point coordinates of testee configuration
  #"dim" to compare the configurations in 2D or the maximum available ("All") dims
  #"pred.dat1" and "pred.dat2" are the predicted responses obtained from the CLPred() function
  #"nsamples" number of samples in original data set
  #"pvar" number of variables original data set
  ###################################################################################
  #Value
  #Returns a data frame with the evaluation measures as stated above.
  ###################################################################################
  ###################################################################################
  
  p <- missbp$p
  n <- missbp$n
  
  tempcomp <- myOPA(missbp, compdat)
  PS <- tempcomp[[1]]
  
  compPred <- CLPpred(CLPs=tempcomp$compCLP, Zs=tempcomp$compZ, p=ncol(compdat), n=nrow(compdat), lvls=tempcomp$complvls, datIN = compdat)
  if(is.data.frame(missbp$X)) 
    {
    datIN <- df2fact(missbp$X)
    } else 
    {
      datIN <- df2fact(missbp$X[[1]])
      }
  GPAPred <-  CLPpred(CLPs=missbp$CLP.GPAbin,Zs=missbp$Z.GPAbin, p=missbp$p, n=missbp$n, lvls=missbp$lvlv[[1]], datIN=datIN)
  
  Z.Target <- tempcomp[[2]]
  CLP.Target <- tempcomp[[3]]
  
  Z.Testee <- missbp$Z.GPAbin
  CLP.Testee <- missbp$CLP.GPAbin
  
  require(stringr)
  
  counter <- 0
  Target <- rbind(Z.Target, CLP.Target)
  Testee <- rbind(Z.Testee, CLP.Testee)
  
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
  
  RMSB <- ((sum(sum((Target-Testee)^2)))/length(Testee))^(0.5)
  AMB <- (sum(sum(abs(Target-Testee))))/length(Testee)

  match.count <- sum(mapply(as.character,compPred)==mapply(as.character,GPAPred))
  
  SP <- match.count/(p*n)

  comp_flat <- apply(compPred,1,str_flatten)
  imp_flat <- apply(GPAPred,1,str_flatten)
  
  RPR <- sum(comp_flat == imp_flat)/nrow(compPred)

  EVALtable <- data.frame(c(PS, SP, RPR, AMB, RMSB))
  colnames(EVALtable)<- c("Evaluation measures")
  rownames(EVALtable)<- c("PS","SP","RPR", "AMB", "RMSB")
  
  missbp$eval <- round(EVALtable,4)
  missbp
}
