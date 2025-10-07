###################################################################################
#' Biplot function
#' 
#' Creates a multiple correspondence analysis (MCA) biplot
#'
#' @param missbp An object of class \code{missbp} obtained from preceding function \code{missmi()}
#' @param Z.col Colour of sample coordinates
#' @param CLP.col Colour of category level point coordinates 
#' @param Z.pch Plotting character of sample coordinates
#' @param CLP.pch Plotting character of category level point coordinates
#' @param Z.cex Size of plotting character for sample points
#' @param CLP.cex Size of plotting character for category level point points 
#' @param title Title of the plot
#'
#' @export
#' 
#' @return
#' \itemize{
#' \item If `compdat = NULL` in \code{\link{evalMeas}}, only a GPAbin biplot will be constructed.
#' \item If a complete data set (`compdat`) was specified in \code{\link{evalMeas}}, two biplots will be constructed: (1) Complete MCA biplot and (2) GPAbin biplot.
#' }
#'  
#' @examples
#' data(implist)
#' missbp <- missmi(implist)|> DRT() |> GPAbin() |> biplFig()
#' 
biplFig <- function (missbp, Z.col="#61223b", CLP.col="#b79962", Z.pch=19, CLP.pch=15, Z.cex=1.5, CLP.cex=1.7, title="") 
{
  CLPs <- missbp$CLP.GPAbin
  Zs <- missbp$Z.GPAbin
  lvls <- missbp$lvls[[1]] #use the first list element for levels, check for other cases (to do)
  
  oldpar <- graphics::par(no.readonly = TRUE) 
  on.exit(graphics::par(oldpar))
  
  grDevices::dev.new()
  
  graphics::par(pty = "s")
  
  #construct two plots if coordinates for complete case is available
  if(is.null(missbp$compCLPs)) {
    plot(rbind(CLPs[,1:2], Zs[,1:2]), pch="", xaxt="n", yaxt="n", xlab="", ylab="", main=title)
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
    } else {
    compCLPs <- missbp$compCLPs
    compZs <- missbp$compZs
    complvls <- missbp$complvls
    
    plot(rbind(compCLPs[,1:2],compZs[,1:2]), pch="", xaxt="n", yaxt="n", xlab="", ylab="", main="Complete biplot")
    graphics::points(compZs, pch=Z.pch, col=Z.col)
    graphics::points(compCLPs, pch=CLP.pch, col=CLP.col)
    
    is.null(complvls)
    {
      graphics::text(compCLPs, cex=0.7, label=rownames(compCLPs), pos=3)
    }
    !is.null(complvls)
    {
      graphics::text(compCLPs, cex=0.7, label=complvls, pos=3)
    }
    grDevices::dev.new()
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
  missbp
}

###################################################################################