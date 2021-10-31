#' Calculate and map matching quality
#'
#'
#'
#' @param x Output from \code{\link{kpoints}} function.
#'
#' @return Plot of the proportion of the study area covered for each value of k,
#' or if only one value of k was used, reports coverage for that solution.
#'
#' @export
#' @examples
#' # Load bioclim data for Target Cells (from rMultivariateMatchingAlgorithms package)
#' data(bioclim)
#' Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(bioclim)

matchingquality <- function(x){

# Standardize variables of interest according to percentage of range
bioclim1 <- cbind(bioclim1a[,1:2],bioclim1a$bioclim_01/0.7,
                  bioclim1a$bioclim_04/42,
                  bioclim1a$bioclim_09/3.3,
                  bioclim1a$bioclim_12/66,
                  bioclim1a$bioclim_15/5.4,
                  bioclim1a$bioclim_18/18.4)
# fix column names
colnames(bioclim1) <- c(colnames(bioclim1a)[1:2],"bioclim_01","bioclim_04","bioclim_09", "bioclim_12","bioclim_15","bioclim_18")

# Calculate distance matrix for all cells:
xdist <- distances(bioclim1[,3:ncol(bioclim1)], id_variable = row.names(bioclim1))
cell_numbers <- rownames(bioclim1)

# Read in subset cells (with cellnumbers as rownames)
centers <- read.csv(paste0(datafolder,"/SubsetCells_siteselection.csv"), row.names = 1)

# Find subset cell that is nearest neighbor to each target cell
neighbors <- nearest_neighbor_search(xdist, k = 1, search_indices = which(rownames(bioclim1) %in% rownames(centers)))

# Create vector of cellnumbers of subset cells matched to each target cell
neighbors2 <- cell_numbers[as.numeric(t(neighbors))]

# Bind to transformed data
bioclim2 <- cbind(bioclim1,neighbors2)

# Use neighbors2 column in bioclim2 to calculate weighted Euclidean distance between target and matched subset cells
d1 <- (bioclim2[as.character(bioclim2[,ncol(bioclim2)]),3]-bioclim2[,3])^2
for (cv in 2:(ncol(bioclim2)-3)){
  d1<- cbind(d1,(bioclim2[as.character(bioclim2[,ncol(bioclim2)]),2+cv]-bioclim2[,2+cv])^2)
}
sum6 <- apply(d1,1, sum)

# Final weighted Euclidean distance between each Target cell and its matched Subset cell (i.e. matching quality variable):
qual <- data.frame(distance = sqrt(sum6))
}
