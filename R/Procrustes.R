#GPA related functions
GPAbin <- function(missbp, G.target=NULL)
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
