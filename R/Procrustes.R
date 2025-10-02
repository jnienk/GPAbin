#' Function to unify coordinates of multiple configurations
#'
#' Combines multiple configurations from dimension reduction solutions applied to multiple imputed data sets
#' 
#' @param missbp An object of class \code{missmi} obtained from preceding function \code{missmi()}
#' @param G.target Target configuration. If not specified the centroid configuration will be used as the target.
#'
#' @return
#' The `missbp` object is appended with the following objects:
#' \item{Z.GPA.list}{List containing the sample coordinates for each MI after GPA}
#' \item{CLP.GPA.list}{List containing the CLPs for each MI after GPA}
#' \item{G.target}{Target configuration}
#' \item{Z.GPAbin}{Sample coordinates for the GPAbin biplot}
#' \item{CLP.GPAbin}{CLPs for the GPAbin biplot}
#' 
#' See also \code{\link{missmi}}, \code{\link{impute}} and \code{\link{DRT}}.
#' 
#' For more detail, refer to Nienkemper-Swanepoel, J., le Roux, N. J., & Gardner-Lubbe, S. (2021). GPAbin: unifying visualizations of multiple imputations for missing values. Communications in Statistics - Simulation and Computation, 52(6), 2666–2685. https://doi.org/10.1080/03610918.2021.1914089.
#' @export
#'
#' @examples
#' data(implist)
#' missbp <- missmi(implist) |> DRT() |> GPAbin()
#'
GPAbin <- function(missbp, G.target=NULL)
{
  m <- missbp$m
  
  z_split <- nrow(missbp$Z[[1]])
  clp_split <- nrow(missbp$CLP[[1]])
  
  coord_set <- vector("list",m)
  for (i in 1:m)
  {
  coord_set[[i]] <- rbind(missbp$Z[[i]], missbp$CLP[[i]])
  }
  tot_rows <- nrow(coord_set[[1]])
  
  #function requires a list
  GPA.out <- GPA(coord_set,G.target=NULL)
  G.target <- GPA.out[[4]]
  Q.list <- GPA.out[[3]]
  s.list <- GPA.out[[2]]
  GPA.list <- GPA.out[[1]]
  
  Z.GPA.list <- vector("list",m)
  CLP.GPA.list <- vector("list",m)
  
  for (i in 1:m)
  {
    Z.GPA.list[[i]] <- GPA.list[[i]][1:z_split,]
    CLP.GPA.list[[i]] <- GPA.list[[i]][((z_split+1):tot_rows),]
  }
  
Z.GPAbin <- Reduce("+",Z.GPA.list)/length(Z.GPA.list)
CLP.GPAbin <- Reduce("+",CLP.GPA.list)/length(CLP.GPA.list)
  
missbp$Z.GPA.list <- Z.GPA.list
missbp$CLP.GPA.list <- CLP.GPA.list
missbp$G.target <- G.target
missbp$Z.GPAbin <- Z.GPAbin
missbp$CLP.GPAbin <- CLP.GPAbin

missbp
}

###################################################################################
###################################################################################
#' Generalised Orthogonal Procrustes Analysis
#' 
#' This function contains the OPA function to compare two configurations and the GPA function for multiple configuration comparisons
#'
#' @param Xk list containing the testee configurations which is updated on #each iteration
#' @param G.target Target configuration. If not specified the centroid configuration will be used as the target
#' @param iter Number of iterations allowed before convergence
#' @param eps Threshold value for convergence of the alogrithm
#'
#' @return
#' \item{Xk.F}{List containing the updated testee configurations}
#' \item{sk.F}{Vector containing the final scaling factors}
#' \item{Qk.F}{List containing the final rotation matrices}
#' \item{Gmat}{Final target configuration}
#' \item{sum.sq}{Final minimised sum of squared distance}
#' 
#' @export
#'
GPA <- function(Xk, G.target=NULL, iter=500, eps=0.001)
{
  OPA.steps <- function(X.mat, Z.mat)
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
    stop(paste("Maximum number of specified iterations reached! Increase iter \n",iter))
  Xk.F <- sapply(1:K, function(k, Xk, sk, Qk) sk[k] * Xk[[k]] %*% Qk[[k]], Xk = Xk.F, sk = sk, Qk = Qk, simplify = F)
  if (is.null(G.target))
  {Gmat <- Xk.F[[1]]
  for(k in 2:K)  Gmat <- Gmat + Xk.F[[k]]
  Gmat <- Gmat/K
  }
  else Gmat <- G.target
  Qk <- sapply(1:K, function(k, Xind, Gmat) OPA.steps(Xind[[k]], Gmat), Xind = Xk.F, Gmat = Gmat, simplify = F)
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

#' Orthogonal Procrustes Analysis
#' 
#' This function performs Orthogonal Procrustes Analysis on centred data
#'
#' @param missbp An object of class \code{missmi} obtained from preceding function \code{missmi()}
#' @param compdat Complete data set, only available for simulated data examples.
#' @param centring Logical argument to apply centering, default is `TRUE`.
#' @param dim Number of dimensions to use in final solutions (`2D` or `All` available dimensions.)
#'
#' @return
#' \item{ProcStat}{Procrustes Statistic}
#' \item{compZ}{Sample coordinates representing the complete data set}
#' \item{compCLP}{Category level point coordinates representing the complete data set}
#' \item{complvls}{Category levels}
#' \item{compdat}{Complete data set, only available for simulated data examples}
#' @export
#' 
OPA <- function(missbp, compdat, centring = TRUE, dim = "2D")
{
  #creating coordinates of complete set
  #DRT on complete data
  temp <- ca::mjca(compdat,lambda="indicator")
  compZ <-temp[[16]]
  compCLP <- temp[[23]]
  complvls <- temp[[6]]
  target <- rbind(compZ, compCLP)
  
  #testee from missbp object
  testee <- rbind(missbp$Z.GPAbin, missbp$CLP.GPAbin)
  
  nCLTar <- nrow(target)
  nCLTes <- nrow(testee)
  
  Tarnam <- complvls#rownames(target)
  Tesnam <- missbp$lvls[[1]]#rownames(testee)

  #finding the CLs that occur in both Target and Testee and deleting the CLs that do #not appear in both in order to obtain a one-to-one comparison
  counter <- 0
  rem <- which(is.na(match(Tarnam,Tesnam)))
  if(is.integer0(rem))
  {
    target <- target
    counter <- counter+1#counts the matched cases
  } else {target <- target[-rem,]}
  
  if(dim=="All")
  {
    pY <- ncol(target)
    pX <- ncol(testee)
    #the maximum number of common columns to use
    colUse <- min(pY,pX)
    target <- target[,1:colUse]
    testee <- testee[,1:colUse]
  } else
    if (dim=="2D")
    {
      target <- target[,1:2]
      testee <- testee[,1:2]
    }
  
  if(!centring)
  {
    testee <- testee
    target <- target
  } else
  {
    testee <- scale(testee,T,F)
    #centre=T, scale=F results are similar to Cox and Cox, Gower and #Dijkersthuis, Borg and Groenen
    target <- scale(target,T,F)
  }
  
  #transformations
  C.mat <- t(target)%*%testee
  svd.C <- svd(C.mat)
  A.mat <- svd.C[[3]]%*%t(svd.C[[2]])
  s.fact <- sum(diag(t(target)%*%testee%*%A.mat))/sum(diag(t(testee)%*%testee))
  #Gower and Dijksterhuis P32
  b.fact <- as.vector(1/nrow(target) * t(target - s.fact * testee %*% A.mat)%*%rep(1,nrow(target)))
  
  X.new <- b.fact + s.fact*testee%*%A.mat
  
  Res.SS <- sum(diag(t(((s.fact*testee%*%A.mat)-target))%*%((s.fact*testee%*%A.mat)-target)))
  PS <- Res.SS/sum(diag(t(target)%*%target))
  RMSB <- ((sum(sum((target-testee)^2)))/length(testee))^(0.5)
  AMB <- (sum(sum(abs(target-testee))))/length(testee)
  
  return(list(ProcStat = PS, RMSB = RMSB, AMB = AMB, compZ=compZ, compCLP=compCLP, complvls=complvls, compdat=compdat))
}
