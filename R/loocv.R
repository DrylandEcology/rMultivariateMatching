#' Leave-one-out cross-validation
#'
#' Use leave-one-out cross-validation of simulated sites to evaluate matching
#' errors. This function takes each Subset cell and finds its nearest neighbor
#' from among the remaining Subset cells using weighted (standardized) Euclidean
#' distance of the matching variables, then calculates the differences between the
#' simulated value of output variables for each Subset cell and the simulated value
#' of output variables from its nearest neighbor.
#'
#' This function can be used for matching achieved with the \code{\link{multivarmatch}}
#' function or matching that uses two-step matching, first with
#' \code{\link{multivarmatch}}, followed by \code{\link{secondaryMatching}}. In
#' case of the latter, `loocv` assumes a very specific experimental design. For
#' each Subsetcell cell, there are five different soil types. There is one site-specific
#' soil type that is unique to each Subset cell and four set soil soil types that
#' were simulated for all Subset cells. In this case, the first step of matching
#' uses only climate variables and the second step of matching identifies the
#' best soil type from among the 5 available for the Subset cell with the best
#' match based on the climate variables.
#'
#'
#' @param matchingvars a data frame that includes all matching variables for the
#' Subset cells. Rownames should correspond to the unique identifiers for each
#' Subset cell. The first two columns correspond to 'x' and 'y' coordinates of
#' the Subset cells (if none exist, use "NA"). The rest of the columns correspond
#' to the matching variables.
#'
#' @param secondaryvars a data frame that includes the secondary matching variables
#' for the Subset cells. The first column should correspond to the unique identifier
#' for each Subset cell, and the next two columns should correspond to the 'x' and
#' 'y' coordinates of the Subset cells (if non exist, use "NA"). The rest of the
#' columns correspond to the secondary matching variables. Only needed if
#' `secondarymatch` is TRUE.
#'
#' @param output_results data frame. Simulation output results for all simulated
#' sites (Subset cells). The first column and the rownames should correspond to
#' the unique identifiers for the Subsetcells. If `secondarymatch` is TRUE,
#' the unique identifiers should be a combination of the unique identifiers used
#' for 'matchingvars' and 'secondaryvars', in the form 'matchingvars_id.treatment'.
#' 'treatment' is a number that designates the treatments, where '1' corresponds
#' to the unique site-specific treatment (i.e., site-specific soils) and subsequent
#' numbers designate other treatments that were applied to all Subset cells
#' (these will correspond to the rownames in `other_treatments`, see below).
#'
#' @param criteria1 single value or vector of length equal to the number of matching variables,
#' where values correspond to the matching criterion for each matching variable
#' in 'matchingvars'. If a single value, this will be used as matching criteria
#' for all variables. Default value is 1, corresponding to using raw data for
#' matching.
#'
#' @param criteria2 single value or vector of length equal to the number of
#' secondary variables, where values correspond to the matching criterion for
#' each secondary variable `secondaryvars`. If a single value, this will be used
#' as matching criteria for all variables. Default value is 1, corresponding to
#' using raw data for matching. Only needed if `secondarymatch` is TRUE.
#'
#' @param secondarymatch boolean. Indicates whether the function should run
#' secondary matching on the Subset cells. Defaults to TRUE
#'
#' @param secondaryvars_id character. Provides the column name for the unique
#' identifiers in the 'secondaryvars' data frame (should be the first column).
#' Defaults to "cellnumbers". Only needed if `secondarymatch` is TRUE.
#'
#' @param reference_treatment character. Designates a number to represent the
#' reference treatment. Default value is '1'. Only needed if `secondarymatch` is
#' TRUE.
#'
#' @param n_neighbors numeric. The number of nearest neighbors to search for among
#' the Subset cells. To achieve leave-one-out cross-validation, this number must be
#' set to 2. The nearest neighbor of each Subset cells is itself, so the second
#' nearest neighbor will correspond the closest non-self neighbor. Default value
#' is 2.
#'
#' @param other_treatments a data frame that gives the secondary variables for the
#' set treatments. The rownames should correspond to unique identifiers for each
#' treatment (e.g., 2-total number of treatments if the `reference_treatment` is '1').
#' Only needed if `secondarymatch` is TRUE.
#'
#' @param ... additional parameters to be passed to \code{\link{multivarmatch}}
#' and/or \code{\link{secondaryMatching}}.
#'
#' @return Data frame of Target cells with coordinates ('x','y'), cellnumber of
#' Target cell ('target_cell'), unique id of matched Subset cell ('subset_cell')
#' and matching quality ('matching_quality'). (If `secondarymatch` is TRUE, it will
#' also include the unique id of the Subset cell matched with secondary matching
#' criteria ('subset_cell_secondary'), and matching quality of this secondary match
#' ('matching_quality_secondary')). Additional columns correspond to the difference
#' between the simulated values of output variables for each Subset cell and the
#' simulated values of output variables from its nearest neighbor. These values
#' can be squared and averaged to get the average squared cross-validated
#' matching error for each output variable.
#'
#'
#' @examples
#' ########################
#' # First, an example where secondarymatch = FALSE
#' # Load targetcells data for Target Cells
#' data(targetcells)
#'
#' # Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#'
#' # Create a mock dataset of output results
#' output_results <- allvars[rownames(subsetcells),c("cellnumbers","bioclim_02",
#'                                                   "bioclim_03","bioclim_16",
#'                                                   "bioclim_17")]
#'
#' # Create dataset of matchingvars for subsetcells
#' subset_matchingvars <- matchingvars[rownames(subsetcells),-1]
#'
#' # Run leave-one-out cross validation of mock output results
#' loocv_results <- loocv(matchingvars = subset_matchingvars,
#'                        output_results = output_results,
#'                        criteria1 = criteria,
#'                        secondarymatch = FALSE, n_neighbors = 2)
#'
#'
#' ########################
#' # Next, an example where secondarymatch = TRUE
#' # Get subsetcells
#' data(subsetcells)
#'
#' # Pull out only matching variables and remove duplicates
#' matchingvars <- subsetcells[,c("site_id","X_WGS84","Y_WGS84","bioclim_01",
#' "bioclim_04","bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
#'
#' # Fix names
#' names(matchingvars) <- c("cellnumbers","x","y",names(matchingvars)[4:9])
#'
#' # Remove duplicates (we will first match on climate only)
#' matchingvars <- matchingvars[!duplicated(matchingvars$cellnumbers),]
#' rownames(matchingvars) <- matchingvars$cellnumbers
#'
#' # Remove cellnumbers column
#' matchingvars <- matchingvars[,-1]
#'
#' # Pull out secondary vars and keep both identifiers
#' secondaryvars <- subsetcells[,c("site_id","X_WGS84","Y_WGS84","sand","clay")]
#'
#' # Fix names
#' names(secondaryvars) <- c("cellnumbers","x","y",names(secondaryvars)[4:5])
#'
#' # Convert sand and clay to percentage from fraction
#' secondaryvars$sand <- secondaryvars$sand*100
#' secondaryvars$clay <- secondaryvars$clay*100
#'
#' # Remove duplicates
#' secondaryvars <- secondaryvars[!duplicated(secondaryvars$cellnumbers),]
#'
#' # Set rownames as cellnumbers
#' rownames(secondaryvars) <- secondaryvars$cellnumbers
#'
#' # Bring in "other" treatments
#' data(setsoiltypes)
#' other_treatments = setsoiltypes
#'
#' # Set original criteria (from first-step matching)
#' criteria1 = c(0.7,42,3.3,66,5.4,18.4)
#'
#' # Calculate criteria for secondary matching
#' criteria2 = c((max(subsetcells$sand,na.rm = T)-
#'              min(subsetcells$sand,na.rm = T))/10*100,
#'              (max(subsetcells$clay,na.rm = T)-
#'              min(subsetcells$clay,na.rm = T))/10*100)
#'
#' # Bring in simulation output results of interest
#' output_results = subsetcells[,c("site_ids","Dryprop","CwetWinter",
#' "CdrySummer","Cwet8","Dryall","Dryany")]
#' rownames(output_results) <- output_results$site_ids
#'
#' # Run leave-one-out cross validation of output results
#' loocv_results <- loocv(matchingvars = matchingvars,
#'                        secondaryvars = secondaryvars,
#'                        output_results = output_results,
#'                        criteria1 = criteria1, criteria2 = criteria2,
#'                        secondarymatch = TRUE,
#'                        secondaryvars_id = "cellnumbers",
#'                        reference_treatment = "1", n_neighbors = 2,
#'                        other_treatments = other_treatments)


loocv <- function(matchingvars, secondaryvars, output_results = NULL,
                  criteria1 = 1, criteria2 = 1, secondarymatch = TRUE,
                  secondaryvars_id = "cellnumbers", reference_treatment = "1",
                   n_neighbors = 2,
                  other_treatments = NULL, ...){

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
                              criteria = criteria2,
                              reference_treatment = "1",
                              raster_template = NULL, saveraster=FALSE,
                              plotraster=FALSE,
                              other_treatments = other_treatments,
                              is_loocv = TRUE, ...)
    # set which column to pull subset cells from
    thissubset = "subset_cell_secondary"
  } else if (!secondarymatch){
    quals2 <- qual1
    # set which column to pull subset cells from
    thissubset = "subset_cell"
  }

  # Now calculate matching errors and add onto matching quality
  for (i in 2:ncol(output_results)){
    quals2[,ncol(quals2)+1] <- output_results[quals2$target_cell, i] - output_results[quals2[,thissubset],i]
    # Fix names
    names(quals2)[ncol(quals2)] <- paste0(names(output_results)[i],"_diff")

  }
  return(quals2)
}

#' Estimate matching errors
#'
#' Use the output from \code{\link{loocv}} to estimate matching errors for each
#' output variable
#'
#'
#' @param loocv_output a data frame produced by the \code{\link{loocv}} function.
#'
#'
#' @return a named vector where items are the estimated matching error for each
#' output variable calculated using leave-one-out cross validation of the
#' Subset cells (simulated sites).
#'
#'
#' @examples
#' ########################
#' # An example where secondarymatch = FALSE
#' # Load targetcells data for Target Cells
#' data(targetcells)
#'
#' # Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#'
#' # Create a mock dataset of output results
#' output_results <- allvars[rownames(subsetcells),c("cellnumbers","bioclim_02",
#'                                                   "bioclim_03","bioclim_16",
#'                                                   "bioclim_17")]
#'
#' # Create dataset of matchingvars for subsetcells
#' subset_matchingvars <- matchingvars[rownames(subsetcells),-1]
#'
#' # Run leave-one-out cross validation of mock output results
#' loocv_results <- loocv(matchingvars = subset_matchingvars,
#'                        output_results = output_results,
#'                        criteria1 = criteria,
#'                        secondarymatch = FALSE, n_neighbors = 2)
#'
#' # Calculate estimates of matching error for output variables
#' estimated_errors <-  cverrors(loocv_output = loocv_results,
#' first_output_column = 6)


cverrors <- function(loocv_output = NULL, first_output_column = NULL){
  results <- vector()
  for (i in first_output_column:ncol(loocv_output)){
  results[i-first_output_column+1] <- sqrt(mean(loocv_output[,i]^2))
  }
  names(results) <- gsub("_diff","",
                         names(loocv_output)[first_output_column:ncol(loocv_output)])
  return(results)
}
