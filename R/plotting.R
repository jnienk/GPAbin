#plotting
###################################################################################
biplFig <- function (missbp, Z.col="#61223b", CLP.col="#b79962", Z.pch=19,
                     CLP.pch=15,Z.cex=1.5,CLP.cex=1.7,title="") 
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
  CLPs <- missbp$CLP.GPAbin
  Zs <- missbp$Z.GPAbin
  lvls <- missbp$lvls[[1]] #use the first list element for levels, check for other cases (to do)
  
  grDevices::dev.new()
  graphics::par(pty="s")
  plot(rbind(CLPs[,1:2],Zs[,1:2]), pch="", xaxt="n", yaxt="n", xlab="", ylab="", main=title)
  graphics::points(Zs, pch=Z.pch, col=Z.col)
  graphics::points(CLPs, pch=CLP.pch, col=CLP.col)
  
  is.null(lvls)
  {
    graphics::text(CLPs, cex=0.7, label=rownames(CLPs), pos=3)
  }
  !is.null(lvls)
  {
    graphics::text(CLPs, cex=0.7, label=lvls, pos=3)
  }
}

###################################################################################