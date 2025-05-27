#Dimension reduction

DRT <- function(missbp, method=c("MCA"))
{
  #in this version only multiple correspondence analysis
  require(ca) 
  
  #if completed data were provided
  if(is.null(missbp$miss_pct)) 
    {
    #continue with imputed data
    data <- missbp$X
  } else{
    data <- missbp$dataimp
    }
  
  m <- missbp$m
  
  Z.list <- vector("list", m)
  CLP.list <- vector("list", m)
  lvls.list <- vector("list", m)
  
  #ensure data is a factor
  for (i in 1:m)
  {
  data[[i]] <- df2fact(data[[i]])
  }
  #remove empty category levels prior to dimension reduction
  data <- rmOneCL(data)
  
  for (i in 1:m)
  {
  temp <- ca::mjca(data[[i]],lambda="indicator")
  Z.list[[i]] <-temp[[16]]
  CLP.list[[i]] <- temp[[23]]
  lvls.list[[i]] <- temp[[6]]
  }
  
  missbp$Z <- Z.list
  missbp$CLP <- CLP.list
  missbp$levels <- lvls.list
  
  missbp
  }