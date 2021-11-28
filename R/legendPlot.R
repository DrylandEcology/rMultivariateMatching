#' Plot a raster with a custom legend below
#'
#' Plots a raster with binned colors and adds a legend below the image.
#'
#'
#' @param x raster. A raster to be plotted.
#'
#' @param thisVariable character. The name of the variable in the raster. Defaults
#' to the name of the raster.
#'
#' @param round_dec numeric. The number of decimal places to round the labels
#' of the legend. Defaults to 0.
#'
#' @param cols vector of colors to use in the legend. Default colors range from
#' yellow to dark brown through 8 distinct colors.
#'
#' @param bks vector of breaks to use in designating which values get assigned to
#' each color in `cols`. Must have one more element than `cols`. Unless `bks` are
#' designated, this vector will be calculated internally from `cols`.
#'
#' @param addpoints boolean. Indicates whether the locations of subset cells
#' should be added to the map as points. Defaults to FALSE.
#'
#' @param matchingQ boolean. Indicates whether the map to be plotted is of
#' matching quality. Defaults to FALSE.
#'
#' @return a plot of the raster with a legend.
#'
#' @author Rachel R. Renne
#'
#' @export
#'
#' @importFrom graphics layout
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom graphics box
#' @importFrom graphics mtext
#' @importFrom graphics polygon
#' @importFrom graphics axis
#' @importFrom graphics legend
#'
#'
#' @examples
#' # Load targetcells data for Target Cells
#' data(targetcells)
#'
#' # Plot the raster with a legend
#' legendPlot(targetcells[[1]])

legendPlot <- function(x = NULL, thisVariable = names(x)[1], round_dec = 0,
                       cols = NULL, bks = NULL, addpoints = FALSE,
                       matchingQ = FALSE){
  if (is.null(x)){
    stop("Verify inputs: 'x' is missing.")
  }
  if (is.null(cols)){
    cols = c("#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#993404","#662506")
  }
  # Determine number of colors:
  n_cols = length(cols)
  if (is.null(bks)){
  bks <- calculateBreaks(x, n_cols = n_cols)
  }
  # Set up axis label locations
  if (!matchingQ){
  axisat <- 0.75
  for (ii in 2:(length(bks))){
    axisat <- c(axisat, (axisat[ii-1]+9.5/length(cols)))
  }
  }

  # Setup layout
  layout.matrix <- matrix(c(1,2),nrow = 2, ncol = 1)
  layout(layout.matrix, widths = c(1,2), heights = c(6,1.5))

  # create figures
  par(mar=c(0,0.5,0.5,0.5), mgp = c(1.5,0.3,0))
  raster::image(x, useRaster = T,
                col = cols,
                breaks = bks,
                xlab = "", ylab ="",
                bty = "n", xaxt = "n",yaxt="n")
    if (addpoints){
    points(subsetcells[,"y"] ~ subsetcells[,"x"], pch = 16, cex = 1, col = "black")
    }
  box()
  par(mar=c(0,0.5,0.2,0.5))
  plot(1:10~1, col = "white", xaxt = "n", yaxt = "n", bty = "n")
  mtext(thisVariable, side = 1, line = -1.01, cex = 1)
  if (!matchingQ){
  for (xx in 1:length(cols)){
    polygon(x = c(axisat[xx],axisat[xx],axisat[xx+1],axisat[xx+1]),
            y = c(7,10,10,7),border = cols[xx], col = cols[xx])
  }
  polygon(x = c(axisat[9],axisat[1],axisat[1],axisat[9]),
          y = c(7,7,10,10), lwd = 1.5)
  par(tcl = -0.3)
  axis(side = 1, pos = 7, at = seq(min(axisat),max(axisat), length.out = length(bks)),
       cex.axis = 0.9, labels = round(bks,round_dec))
  } else if (matchingQ){
    legend("top" , legend = c("0 to 0.5","0.5 to 1","1 to 1.5",paste0("1.5 to ",ceiling(max(raster::values(x), na.rm=TRUE)))), fill = cols,
           bty = "n", cex = 1, ncol = 4, x.intersp = 0.5)
  }
}


#' A helper function for \code{\link{legendPlot}}
#'
#' Creates equal sized bins for the legend and display of the raster.
#'
#'
#' @param x raster.
#'
#' @param n_cols the number of colors for the legend. Takes input from
#' \code{\link{legendPlot}}
#'
#' @return a vector of breaks for the legend.
#'
#' @importFrom stats na.omit
#'
#' @author Rachel R. Renne
#'
#' @export

calculateBreaks <- function(x, n_cols){
  minx = min(na.omit(raster::values(x)))
  maxx = max(na.omit(raster::values(x)))
  rangex = maxx-minx
  bks = minx
  for (cx in 2:(n_cols)){
    bks[cx] <- bks[cx-1]+rangex/n_cols
  }
  bks[length(bks)+1] = maxx
  return(bks)
}
