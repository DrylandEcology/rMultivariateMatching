#' Estimate matching error
#'
#' Use leave-one-out cross-validation of simulated sites to evaluate matching
#' errors
#'
#'
#' @param
#'
#' @param
#'
#' @param
#'
#' @param
#'
#' @param
#'
#' @return
#'
#'
#' @examples
#' # Get subsetcells
#' data(subsetcells)
#' # Pull out only matching variables and remove duplicates
#' matchingvars <- subsetcells[,c("site_id","X_WGS84","Y_WGS84","bioclim_01","bioclim_04",
#'                                "bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
#'
#' # Fix names
#' names(matchingvars) <- c("cellnumbers","x","y",names(matchingvars)[4:9])
#' # Remove duplicates (we will first match on climate only)
#' matchingvars <- matchingvars[!duplicated(matchingvars$cellnumbers),]
#' rownames(matchingvars) <- matchingvars$cellnumbers
#' # Remove cellnumbers column
#' matchingvars <- matchingvars[,-1]
#'
#' # Pull out secondary vars and keep both identifiers
#' secondaryvars <- subsetcells[,c("site_id","X_WGS84","Y_WGS84","sand","clay")]
#' # Fix names
#' names(secondaryvars) <- c("cellnumbers","x","y",names(secondaryvars)[4:5])
#' # Convert sand and clay to percentage from fraction
#' secondaryvars$sand <- secondaryvars$sand*100
#' secondaryvars$clay <- secondaryvars$clay*100
#' # Remove duplicates
#' secondaryvars <- secondaryvars[!duplicated(secondaryvars$cellnumbers),]
#' rownames(secondaryvars) <- secondaryvars$cellnumbers
#'
#' # Bring in "other" treatments
#' data(setsoiltypes)
#' other_treatments = setsoiltypes
#'
#' # Set original criteria
#' criteria1 = c(0.7,42,3.3,66,5.4,18.4)
#' # Calculate criteria for secondary matching
#' criteria2 = c((max(subsetcells$sand,na.rm = T)-min(subsetcells$sand,na.rm = T))/10*100,
#'              (max(subsetcells$clay,na.rm = T)-min(subsetcells$clay,na.rm = T))/10*100)
#'
#' # Bring in simulation output results of interest
#' output_results = subsetcells[,c("site_ids","Dryprop","CwetWinter","CdrySummer",
#'                                 "Cwet8","Dryall","Dryany")]
#' rownames(output_results) <- output_results$site_ids


loocv <- function(matchingvars, secondaryvars, ouput_results = NULL,
                  criteria1 = 1, criteria2 = 1, secondarymatch = TRUE,
                  secondaryvars_id = "cellnumbers", reference_treatment = "1",
                   n_neighbors = 2,
                  other_treatments = NULL,
                  raster_template = NULL,subset_in_target = FALSE,
                  saveraster=FALSE,plotraster=FALSE, ...){

  # Standardize variables of interest
  stdvars <- matchingvars[,c("x","y")]
  for (i in which(colnames(matchingvars) != c("x","y"))){
  stdvars <- cbind(stdvars, matchingvars[,i]/criteria1[i-2])
  }

  # fix colnames
  colnames(stdvars) <- c("x","y",colnames(matchingvars)[3:ncol(matchingvars)])

  # Calculate distance matrix for all cells:
  xdist <- distances::distances(stdvars[,3:ncol(stdvars)], id_variable = row.names(stdvars))
  cell_numbers <- rownames(stdvars)

  # Find nearest 2 neighbor (1st will be self, second will be nearest non-self cell)
  neighbors <- distances::nearest_neighbor_search(xdist, k = n_neighbors, search_indices = c(1:nrow(matchingvars)))

  # Create vector of cellnumbers of subset cells matched to each target cell
  neighbors2 <- cell_numbers[as.numeric(t(neighbors[nrow(neighbors),]))]

  # Bind to transformed data
  stdvars2 <- cbind(stdvars,as.numeric(t(neighbors[nrow(neighbors),])))

  # Use neighbors2 column in stdvars2 to calculate weighted squared Euclidean
  # distance between target and matched subset cells
  d1 <- (stdvars2[stdvars2[,ncol(stdvars2)],3]-stdvars2[,3])^2
  for (cv in 2:(ncol(stdvars2)-3)){
  d1<- cbind(d1,(stdvars2[stdvars2[,ncol(stdvars2)],2+cv]-stdvars2[,2+cv])^2)
  }
  sum6 <- apply(d1,1, sum)

  # Final weighted Euclidean distance between each Target cell and its matched
  # Subset cell (i.e. matching quality variable):
  qual1 <- data.frame(x = stdvars2[,"x"], y = stdvars2[,"y"],
                    target_cell = rownames(stdvars2),
                    subset_cell = neighbors2, matching_quality = sqrt(sum6))
  rownames(qual1) <- rownames(stdvars)

  if (secondarymatch){
  # Run secondary matching on soils data
  quals2 <- secondaryMatching(secondaryvars = secondaryvars, matches = qual1,
                              subsetcells = secondaryvars,
                              subsetcells_id = secondaryvars_id,
                              subset_in_target = FALSE, criteria = criteria2,
                              reference_treatment = "1",
                              raster_template = raster_template,
                              other_treatments = other_treatments,
                              is_loocv = TRUE, ...)


  } else if (!secondarymatch){
    quals2 <- qual1
  }

  # Now calculate matching errors and add onto matching quality
  for (i in 2:ncol(output_results)){
    quals2[,ncol(quals2)+1] <- output_results[quals2$target_cell, i] - output_results[quals2$subset_cell_secondary,i]
  }
  # Fix names
  names(quals2)[8:ncol(quals2)] <- paste0(names(output_results)[2:ncol(output_results)],"_diff")
  return(quals2)
}
