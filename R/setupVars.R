#' Set up input data for Target cells
#'
#' Prepare data frame of potential matching variables for Target cells from raster files
#'
#'
#' @param x RasterStack that covers the area of interest and includes all potential matching variables.
#' Large datasets may be difficult to handle if memory is limited and it may work better to run the function
#' for several RasterStacks (removing previous RasterStacks before reading in the next one)
#' and combine the output into one final data frame. Note: rasters must all have the same extent,
#' resolution, and crs.
#'
#' @return data frame with a 'cellnumbers' column, x and y coordinates, and values for each of the
#' rasters in the RasterStack
#'
#' @export
#' @examples
#' # Load bioclim data for Target Cells (from rMultivariateMatchingAlgorithms package)
#' data(bioclim)
#' # Create data frame of potential matching variables for Target Cells
#' y <- makeInputdata(bioclim)

# Function to format rasters of potential matching variables for Target Cells
makeInputdata <- function(x){
  # Convert rasters to points
  x1 <- raster::rasterToPoints(x)
  # Get cellnumbers:
  cellnumbers  <- raster::extract(x[[1]], y=x1[,1:2], cellnumbers = T, sp = T)
  # Bind cellnumbers (first column only) to rest of variables (and coordinates)
  x2 <- cbind(cellnumbers[,1],x1)
  # Fix column names
  colnames(x2)[1] <- c("cellnumbers")
  # Make cellnumbers the rownames also
  row.names(x2) <- x2[,1]
  as.data.frame(x2)
}





#' Select matching variables
#'
#' Create a new factor from two existing factors, where the new factor's levels
#' are the union of the levels of the input factors.
#'
#' @param x data frame created using \code{\link{makeInputdata}} or formatted similarly.

#'
#' @return factor
#' @export
#' @examples
#' # Load bioclim data for Target Cells (from rMultivariateMatchingAlgorithms package)
#' data(bioclim)
#' # Create data frame of potential matching variables for Target Cells
#' y <- makeInputdata(bioclim)

data(bioclim)
# Create data frame of potential matching variables for Target Cells
x <- bioclim[[1]]

areas <- area(wydrylands) # Gives approx. area in km2
rm(wydrylands)

# Read in bioclim and soils data
bioclim <- readRDS(paste0(datafolder,"/TargetCells_Bioclim+soils_WY.csv"))
# Note: cellnumbers in this dataframe correspond to the unique cellnumbers associated with each
# target cell. These serve as unique identifiers throughout the code and can
# be extracted from the raster datasets that define the target cells using the
# extract function in the raster package with "cellnumbers" = TRUE

# Find total area of study area by extracting areas using cellnumbers
totalarea <- round(sum(extract(areas, as.numeric(row.names(bioclim)))))
# 246771 km2

# limit to just the 8 variables of interest
# Selected to capture the major drivers of plant community structure an ecohydrology in drylands (from Renne et al. 2021):
# "bioclim_01" = mean annual temperature
# "bioclim_04" = temperature seasonality
# "bioclim_09" = mean temperature of the driest quarter
# "bioclim_12" = mean annual precipitation
# "bioclim_15" = precipitation seasonality
# "bioclim_18" = precipitation of the warmest quarter
# "clay" = depth-weighted percentage clay
# "sand" = depth-weighted percentage sand
bioclim1a <- data.frame(bioclim[,c(2,3,4,7,12,15,18,21,23,24)])

#################################################################
# Step 2: Examine matching variables to help determine matching critera:
summary(bioclim1a)

# Look at 5th and 95th percentiles
apply(bioclim1a[,3:10], 2, FUN = function(x){quantile(x, probs = c(0.05,.95), na.rm = T)})

# Find range of each
apply(bioclim1a[,3:10],2,function(x){(max(x)-min(x))})

# Compare 5% of range to 10% of range for each variable
apply(bioclim1a[,3:10],2,function(x){(max(x)-min(x))*0.05})
apply(bioclim1a[,3:10],2,function(x){(max(x)-min(x))*0.1})



#' Determine matching criteria
#'
#' Create a new factor from two existing factors, where the new factor's levels
#' are the union of the levels of the input factors.
#'
#' @param a factor
#' @param y any raster layer of the study area, with extent, resolution, and crs matching
#' rasters used to create x.
#'
#' @return factor
#' @export
#' @examples
#' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])
