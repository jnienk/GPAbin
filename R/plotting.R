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
    graphics::par(pty = "s")
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


#' ggBiplot function
#' 
#' Creates a multiple correspondence analysis (MCA) biplot in `ggplot`.
#'
#' @param missbp An object of class \code{missbp} obtained from preceding function \code{missmi()}
#' @param Z.col Colour of sample coordinates
#' @param CLP.col Colour of category level point coordinates 
#' @param Z.pch Plotting character of sample coordinates
#' @param CLP.pch Plotting character of category level point coordinates
#'
#' @returns
#' \itemize{
#' \item {plot}{Returns a GPAbin ggplot biplot.}
#' \item {plotC}{Returns an MCA ggplot biplot of the complete data set, if provided.}
#' }
#' @export
#'
#' @examples
#' data(implist)
#' data(compdat)
#' missbp <- missmi(implist)|> DRT() |> GPAbin() |> 
#' evalMeas(compdat = compdat) |> ggbiplFig()
#' 
#' ### GPAbin biplot
#' 
#' missbp$plot
#' 
#' ### MCA biplot
#' 
#' missbp$plotC
#' 
ggbiplFig <- function (missbp, Z.col="#61223b", CLP.col="#b79962", Z.pch=19, CLP.pch=15) 
{
  Zs <- as.data.frame(missbp$Z.GPAbin)[,1:2]
  lvls <- missbp$lvls[[1]] #use the first list element for levels, check for other cases (to do)
  CLPs <- as.data.frame(missbp$CLP.GPAbin) |> dplyr::mutate(source = lvls)
  
  all_coords <- dplyr::bind_rows(Zs |> dplyr::select("V1","V2"), 
                          CLPs |> dplyr::select("V1","V2", "source"))
  
  #construct two plots if coordinates for complete case is available
    all_coords <- dplyr::bind_rows(Zs |> dplyr::select("V1","V2"), 
                                   CLPs |> dplyr::select("V1","V2", "source"))
    
    missbp$plot <- ggplot2::ggplot() +
      ggplot2::geom_blank(data = all_coords, ggplot2::aes(x = V1, y = V2)) +
      ggplot2::geom_point(data = Zs, ggplot2::aes(x = V1, y = V2)
                 , size = 2.5, colour =  Z.col, shape = Z.pch) +
      ggplot2::geom_point(data = CLPs, ggplot2::aes(x = V1, y = V2),
                 colour = CLP.col, size = 3, shape = CLP.pch) +
      ggrepel::geom_text_repel(data = CLPs,
                               ggplot2::aes(V1, V2, label = source),
                      box.padding = 0.6,
                      point.padding = 0.4,
                      force = 1.5) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(title = "GPAbin biplot", x = NULL, y = NULL) +
      ggplot2::theme_void()
    
    if(!is.null(missbp$compCLPs)) {
    compZs <- as.data.frame(missbp$compZs)[,1:2]
    complvls <- missbp$complvls
    compCLPs <- as.data.frame(missbp$compCLPs) |> dplyr::mutate(source = complvls)
    
    all_coordsC <- dplyr::bind_rows(compZs |> dplyr::select("V1","V2"), 
                            compCLPs |> dplyr::select("V1","V2", "source"))
    
    missbp$plotC <- ggplot2::ggplot() +
      ggplot2::geom_blank(data = all_coordsC, ggplot2::aes(x = V1, y = V2)) +
      ggplot2::geom_point(data = compZs, ggplot2::aes(x = V1, y = V2)
                 , size = 2.5, colour = Z.col, shape = Z.pch) +
      ggplot2::geom_point(data = CLPs, ggplot2::aes(x = V1, y = V2),
                 colour = CLP.col, size = 3, shape = CLP.pch) +
      ggrepel::geom_text_repel(data = compCLPs,
                               ggplot2::aes(V1, V2, label = source),
                      box.padding = 0.6,
                      point.padding = 0.4,
                      force = 1.5) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(title = "MCA biplot", x = NULL, y = NULL) +
      ggplot2::theme_void()
  }  
  missbp
}