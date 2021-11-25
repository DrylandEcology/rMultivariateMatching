#' Identify matches, calculate and map matching quality
#'
#'
#' Identifies matches from Subset cells for all Target cells, then calculates
#' matching quality (weighted Euclidean distance between Target and matched
#' Subset cells), with options to plot map of matching quality and save a raster
#' of matching quailty.
#'
#' @param matchingvars data frame generated using \code{\link{makeInputdata}} or
#' formatted such that: column 1 and rownames are 'cellnumbers' extracted using the
#' \code{\link[raster-extract]{raster::extract()}} function, columns 2 and 3 correspond to x and y
#' coordinates, and additional columns correspond to potential matching variables
#' extracted using the \code{\link[raster-rasterToPoints]{raster::rasterToPoints()}} function. These data
#' represent Target cells.
#'
#' @param subsetcells if `subset_in_target` is TRUE, this should be a data frame
#' of coordinates (expects coordinates in columns named 'x' and 'y') for Subset
#' cells. May be extracted from output from \code{\link{kpoints}} function or
#' provided separately. Row names should be unique identifiers for each point
#' (unique means no repeats in rownames of subsetcells if `subset_in_target` is
#' TRUE). If `subset_in_target` is FALSE, this should be a data frame of subset
#' cells with column names corresponding exactly to those in `matchingvars` and
#' row names should be unique identifiers (unique means no repeats among all
#' row names in targetcells and matchingvars if `subset_in_target` is FALSE).
#' See `subset_in_target`.
#'
#' @param matchingvars_id character or numeric. Refers to the column in
#' `matchingvars`that provides the unique identifiers for Target cells. Defaults
#' to "cellnumbers", which is the unique ID column created by \code{\link{makeInputdata}}.
#'
#' @param subsetcells_id character or numeric, but must be composed of numbers
#' and convertable to numeric. Refers to the column in `subsetcells`that provides
#' the unique identifiers for Subset cells. When `subset_in_target` is TRUE,
#' these ids must be unique from `matchingvars_ids`. Note that if there are
#' repeats between the`matchingvars_id`s and the `subsetcells_id`s, you can paste
#' "00" before the `subsetcells_id`s to ensure they are unique from the
#' `matchingvars_id`s. Defaults to NULL.
#'
#' @param criteria single value or vector of length equal to the number of matching variables,
#' where values corresponds to the matching criterion for each matching variable
#' in x. If a single value, this will be used as matching criteria for all variables.
#' Default value is 1, corresponding to using raw data for matching.
#'
#' @param n_neighbors numeric. The number of neighbors to search for in matching.
#' Default value is 1 and this setting should be used for matching. Option for 2+
#' is only included for leave-one-out cross-validation \code{\link{loocv}}.
#'
#' @param addpoints boolean. Indicates if Subset cells should be added to the plot
#' as points. Defaults to FALSE.
#'
#' @param raster_template one of the raster layers used for input data.
#'
#' @param subset_in_target boolean. Indicates if Subset cells have been selected
#' from Target cells using \code{\link{kpoints}} function.
#'
#' @param saveraster boolean. Indicates if raster of matching quality should be
#' saved to file.
#'
#' @param filepath provides path for location where raster will be saved. Defaults
#' to working directory.
#'
#' @param plotraster boolean. Indicates if raster should be plotted to a map.
#'
#' @return Data frame of Target cells with coordinates ('x','y'), cellnumber of
#' Target cell ('target_cell'), unique id of matched Subset cell ('subset_cell')
#' and matching quality ('matching_quality'). Will save a raster of matching
#' quality if `saveraster` is TRUE and plot a map of matching quality if
#' `plotraster` is TRUE.
#'
#' @examples
#' # Load targetcells data for Target Cells
#' data(targetcells)
#'
#' Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#'
#' # Restrict data to matching variables of interest
#' matchingvars <- allvars[,c("cellnumbers","x","y","bioclim_01","bioclim_04",
#' "bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
#'
#' # Create vector of matching criteria
#' criteria <- c(0.7,42,3.3,66,5.4,18.4)
#'
#' # Find solution for k = 200
#' # Note: n_starts should be >= 10, it is 1 here to reduce run time.
#' results1 <- kpoints(matchingvars,criteria = criteria,klist = 200,
#' n_starts = 1,min_area = 50,iter = 50,raster_template = targetcells[[1]])
#'
#' ###################################
#' # First an example where subset_in_target = TRUE
#' # Get points from solution to kpoints algorithm
#' subsetcells <- results1$solutions[[1]]
#'
#' # Find matches and calculate matching quality
#' quals <- multivarmatch(matchingvars, subsetcells, criteria = criteria,
#'                         matchingvars_id = "cellnumbers", addpoints = FALSE,
#'                         raster_template = targetcells[[1]],
#'                         subset_in_target = TRUE)
#'
#' ###################################
#' # Now an example where subset_in_target is FALSE
#' # Get Subset cells data
#' data(subsetcells)
#'
#' # Remove duplicates (representing cells with same climate but different soils--
#' # we want to match on climate only)
#' subsetcells <- subsetcells[!duplicated(subsetcells$site_id),]
#'
#' # Pull out matching variables only, with site_id that identifies unique climate
#' subsetcells <- subsetcells[,c("site_id","X_WGS84","Y_WGS84","bioclim_01",
#' "bioclim_04","bioclim_09","bioclim_12",
#' "bioclim_15","bioclim_18")]
#'
#' # Ensure that site_id will be values unique to subsetcells
#' subsetcells$site_id <- paste0("00",subsetcells$site_id)
#'
#' # Find matches and calculate matching quality
#' quals <- multivarmatch(matchingvars, subsetcells, criteria = criteria,
#'                          matchingvars_id = "cellnumbers",
#'                          subsetcells_id = "site_id",
#'                          raster_template = targetcells[[1]],
#'                          subset_in_target = FALSE, addpoints = FALSE)


multivarmatch <- function(matchingvars,subsetcells,matchingvars_id = "cellnumbers",
                            subsetcells_id = NULL, criteria = 1, n_neighbors = 1,
                            raster_template,subset_in_target = TRUE,
                            saveraster=FALSE,plotraster=TRUE,filepath=getwd(),
                            overwrite = FALSE, ...){
  # If no standardization or if standardization is the same for all matching variables
  if (length(criteria)==1){
    criteria = rep(criteria, (ncol(matchingvars)-3))
  }
  # If incorrect number of matching criteria provided, stop fuunction.
  if (length(criteria) > 1 && length(criteria) != (length(names(matchingvars))-3)){
    stop("Number of matching criteria unequal to number of matching variables; cannot complete matching.")
  }

  # verify matchingvars data frame has proper setup:
  if (is.null(matchingvars$x) | is.null(matchingvars$y)) {
    stop("Coordinates not found as columns 'x' and 'y' in matchingvars data frame.")
  }
  if (max(which(colnames(matchingvars) %in% c('x','y'))) > 3){
    stop("Coordinates not found in correct columns of matchingvars.")
  }
  if (which(colnames(matchingvars) == matchingvars_id) > 3 | is.null(matchingvars_id)){
    stop("matchingvar_id not found in correct column of matchingvars.")
  }
  for (ci in 4:ncol(matchingvars)){
    if (!is.numeric(matchingvars[,ci])){
      stop("Matching variable '",colnames(matchingvars)[ci],"' is not numeric.")
    }
  }

  # If subset_in_target is FALSE (subset cells provided in separate data frame)
  if (!subset_in_target){
    if (is.null(subsetcells_id)){
      stop("Missing 'subsetcells_id', cannot complete matching.")
    }
    # Coerce id columns to character for comparison:
    matchingvars[,matchingvars_id] <- as.character(matchingvars[,matchingvars_id])
    subsetcells[,subsetcells_id] <- as.character(subsetcells[,subsetcells_id])

    # Check that id variable columns are all unique
    if (sum(matchingvars[,matchingvars_id] %in% subsetcells[,subsetcells_id]) > 0){
      warning("Non-unique identifiers in subsetcells and matchingvars, matching may fail.")
    }

    # Set id columns to rownames
    rownames(matchingvars) <- matchingvars[,matchingvars_id]
    rownames(subsetcells) <- subsetcells[,subsetcells_id]

    matchingvars1 <- rbind(as.matrix(matchingvars[,-which(colnames(matchingvars) == matchingvars_id)]),
                           as.matrix(subsetcells[,-which(colnames(subsetcells) == subsetcells_id)]))

  } else if (subset_in_target) { # If subset_in_target is TRUE (subset cells are a subset of target cells)
    matchingvars1 <- matchingvars[,-which(colnames(matchingvars) == matchingvars_id)]
  }

  # Standardize variables of interest
  stdvars <- matchingvars1[,c("x","y")]
  for (i in which(colnames(matchingvars1) != c("x","y"))){
    stdvars <- cbind(stdvars, matchingvars1[,i]/criteria[i-2])
  }

  # fix colnames
  colnames(stdvars) <- c("x","y",colnames(matchingvars1)[3:ncol(matchingvars1)])

  # Calculate distance matrix for all cells:
  xdist <- distances::distances(stdvars[,3:ncol(stdvars)], id_variable = row.names(stdvars))
  cell_numbers <- rownames(stdvars)

  # Find subset cell that is nearest neighbor to each target cell
  if (subset_in_target){
  neighbors <- distances::nearest_neighbor_search(xdist, k = n_neighbors, search_indices = which(rownames(matchingvars1) %in% rownames(subsetcells)))
  }else if (!subset_in_target){
  neighbors <- distances::nearest_neighbor_search(xdist, k = n_neighbors, search_indices = (nrow(matchingvars)+1):nrow(matchingvars1))
  }

  # Create vector of cellnumbers of subset cells matched to each target cell
  neighbors2 <- cell_numbers[as.numeric(t(neighbors[nrow(neighbors),]))]

  # Bind to transformed data
  stdvars2 <- cbind(stdvars,as.numeric(t(neighbors[nrow(neighbors),])))

  # Use neighbors2 column in stdvars2 to calculate weighted Euclidean distance between target and matched subset cells
  d1 <- (stdvars2[stdvars2[,ncol(stdvars2)],3]-stdvars2[,3])^2
  for (cv in 2:(ncol(stdvars2)-3)){
    d1<- cbind(d1,(stdvars2[stdvars2[,ncol(stdvars2)],2+cv]-stdvars2[,2+cv])^2)
  }
  sum6 <- apply(d1,1, sum)

  # Final weighted Euclidean distance between each Target cell and its matched Subset cell (i.e. matching quality variable):
  qual <- data.frame(x = stdvars2[,"x"], y = stdvars2[,"y"], target_cell = rownames(stdvars2),
                     subset_cell = neighbors2, matching_quality = sqrt(sum6))
  rownames(qual) <- rownames(stdvars)

  # Change subset_cells back to numeric
  qual$subset_cell <- as.character(as.numeric(qual$subset_cell))

  # Limit qual to just target cells (we don't need subset cells included if subset_in_target is FALSE)
  qual <- qual[1:nrow(matchingvars),]

  if (saveraster | plotraster){
  # Create spatial points dataframe
  ptsx <- sp::SpatialPointsDataFrame(qual[,1:2], data = qual,
                                     proj4string = raster::crs(raster_template))

  # Create raster of qual (matching quality) using wydry as a template
  r <- raster::rasterize(ptsx, raster_template, field = qual$matching_quality, fun = mean)

  if (saveraster){
    raster::writeRaster(r,paste0(filepath,"/Matchingquality.tif"),
                        overwrite = overwrite)
  }

  if (plotraster){
  # Designate colors, breaks, and plotting settings
  cols = rev(c("#d7191c","#fdae61","#abd9e9","#2c7bb6"))
  bks = c(0,0.5,1,1.5,5)
  thisVariable = "Matching quality"
  legendPlot(r, bks = bks, cols = cols, thisVariable = "Matching quality",matchingQ = TRUE, ...)
  }
  }
  return(qual)
}




#' Evaluate matching for additional variables
#'
#' Calculate the standard deviation of differences between Subset and Target cells
#' for a set of variables relevant to the project. This is most informative if it
#' includes variables not used for matching.
#'
#' @param allvars data frame generated using \code{\link{makeInputdata}} or
#' formatted such that: column 1 and rownames are 'cellnumbers' extracted using the
#' \code{\link[raster-extract]{raster::extract()}} function, columns 2 and 3
#' correspond to x and y coordinates, and additional columns correspond to various
#' variables (which can include matching variables) that have been extracted to
#' points using the \code{\link[raster-rasterToPoints]{raster::rasterToPoints()}}
#' function. These data represent Target cells (and may also represent Subset
#' cells if `subset_in_target` is TRUE).
#'
#' @param subsetcells if `subset_in_target` is TRUE, this should be a data frame
#' of coordinates (expects coordinates in columns named 'x' and 'y') for Subset
#' cells. May be extracted from output from \code{\link{kpoints}} function or
#' provided separately. Row names should be unique identifiers for each point
#' (unique means no repeats in rownames of subsetcells if `subset_in_target` is
#' TRUE). If `subset_in_target` is FALSE, this should be a data frame of subset
#' cells with column names corresponding exactly to those in `matchingvars` and
#' row names should be unique identifiers (unique means no repeats among all
#' row names in targetcells and matchingvars if `subset_in_target` is FALSE).
#' See `subset_in_target`.
#'
#' @param matchingvars_id character or numeric. Refers to the column in
#' `matchingvars`that provides the unique identifiers for target cells. Defaults
#' to "cellnumbers", which is the unique ID column created by \code{\link{makeInputdata}}.
#'
#' @param subsetcells_id character or numeric, but must be composed of numbers
#' and convertable to numeric. Refers to the column in `subsetcells`that provides
#' the unique identifiers for Subset cells. When `subset_in_target` is TRUE,
#' these ids must be unique from `matchingvars_ids`. Note that if there are
#' repeats between the`matchingvars_id`s and the `subsetcells_id`s, you can paste
#' "00" before the `subsetcells_id`s to ensure they are unique from the
#' `matchingvars_id`s. Defaults to NULL.
#'
#' @param matches data frame output from the \code{\link{multivarmatch}} or
#' \code{\link{secondaryMatching}} functions.
#'
#' @param secondarymatch boolean. Indicates if the `matches` data frame comes
#' from the \code{\link{secondaryMatching}} function.
#'
#' @param quality_name character. Name of the column in the `matches` data frame
#' that contains the matching quality variable to use to evaluate matching
#' "matching_quality" or "matching_quality_secondary".
#' Defaults to "matching_quality"
#'
#' @param subset_in_target boolean. Indicates if Subset cells have been selected
#' from Target cells using \code{\link{kpoints}} function
#'
#' @param exclude_poor_matches boolean. Indicates if Target cells with poor
#' matching quality (large weighted Euclidean distance to matched Subset cell)
#' should be excluded from calculations. Defaults to TRUE.
#'
#' @param matching_distance numeric. Gives the maximum allowable matching quality
#' value (weighted Euclidean distance) between Target and Subset cells. Default
#' value is 1.5.
#'
#' @param plot_diffs boolean. Indicates whether a barplot of differences should
#' be displayed.
#'
#' @return Data frame of the standard deviation of differences between Target and
#' their matched Subset cells for all variables supplied in `allvars` data frame.
#' The first row corresponds to the standard deviation of differences between
#' Target and Subset cells for all cells and the second row corresponds to the
#' standard deviation of differences between Target and Subset cells for only those
#' Target cells with matching quality <= `matching_distance`. Units are the same
#' as the units for each variable in `allvars`.
#'
#' @examples
#' # Load targetcells data for Target Cells (from rMultivariateMatchingAlgorithms package)
#' data(targetcells)
#'
#' # Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#'
#' # Subset to include only matching variables
#' matchingvars <- allvars[,c("cellnumbers","x","y","bioclim_01","bioclim_04",
#' "bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
#'
#' # Create vector of matching criteria
#' criteria <- c(0.7,42,3.3,66,5.4,18.4)
#'
#' # Find solution for k = 200
#' # Note: n_starts should be >= 10, it is 1 here to reduce run time.
#' results1 <- kpoints(matchingvars,criteria = criteria,
#' klist = 200,n_starts = 1,min_area = 50,iter = 50,
#' raster_template = targetcells[[1]])
#'
#'
#' ###################################
#' # First an example where subset_in_target = TRUE
#' # Get points from solution to kpoints algorithm
#' subsetcells <- results1$solutions[[1]]
#'
#' # Find matches and calculate matching quality
#' quals <- multivarmatch(matchingvars, subsetcells,
#' matchingvars_id = "cellnumbers", raster_template = targetcells[[1]],
#' subset_in_target = TRUE)
#'
#' # Run evaluateMatching
#' sddiffs <- evaluateMatching(allvars = allvars, matches = quals,
#'                             matchingvars_id = "cellnumbers",
#'                             secondarymatch = FALSE,
#'                             subset_in_target = TRUE, matching_distance = 1.5,
#'                             plot_diffs = TRUE)
#'
#' ###################################
#' # Now an example where subset_in_target is FALSE
#' # Get points from solution to kpoints algorithm
#' data(subsetcells)
#'
#' # Remove duplicates (representing cells with same climate but different
#' soils--we want to match on climate only)
#' subsetcells <- subsetcells[!duplicated(subsetcells$site_id),]
#'
#' # Pull out matching variables only, with site_id that identifies unique climate
#' subsetcells1 <- subsetcells[,c("site_id","X_WGS84","Y_WGS84","bioclim_01",
#' "bioclim_04","bioclim_09","bioclim_12",
#' "bioclim_15","bioclim_18")]
#'
#' # Ensure that site_id will be values unique to subsetcells
#' subsetcells1$site_id <- paste0("00",subsetcells$site_id)
#'
#' # Find matches and calculate matching quality
#' quals <- multivarmatch(matchingvars, subsetcells=subsetcells1,
#' matchingvars_id = "cellnumbers",subsetcells_id = "site_id",
#' raster_templat = targetcells[[1]], subset_in_target = FALSE)
#'
#' # Get all variables for Subset cells now:
#' subsetcells <- subsetcells[,c("site_id","X_WGS84","Y_WGS84",
#' names(allvars)[4:22])]
#'
#' # Run evaluateMatching
#' sddiffs <- evaluateMatching(allvars = allvars[,c(1:22)],
#'                             subsetcells = subsetcells,
#'                             secondarymatch = TRUE,
#'                             quality_name = "matching_quality",
#'                             matches = quals,
#'                             matchingvars_id = "cellnumbers",
#'                             subsetcells_id = "site_id",
#'                             subset_in_target = FALSE,
#'                             matching_distance = 1.5,
#'                             plot_diffs = TRUE)

evaluateMatching <- function(allvars = NULL, subsetcells = NULL,
                             matches = NULL, secondarymatch = FALSE,
                             quality_name = "matching_quality",
                             matchingvars_id = "cellnumbers",
                             subsetcells_id = NULL,subset_in_target = TRUE,
                             matching_distance = 1.5, plot_diffs = TRUE){
  if (is.null(allvars) | is.null(matches)){
    stop("Verify inputs: 'allvars' or 'matches' is missing.")
  }
  if (sum(names(allvars)[1:3] %in% c('x','y','cellnumbers')) < 3){
    stop("Verify format of the first three columns in 'allvars'. See documentation for details.")
  }
  for (ci in 4:ncol(matchingvars)){
    if (!is.numeric(matchingvars[,ci])){
      stop("Variable '",colnames(allvars)[ci],"' is not numeric.")
    }
  }
  # Modify allvars and matches to exclude missing data
  allvars <- allvars[complete.cases(allvars),]
  matches <- matches[complete.cases(matches),]
  if (sum(rownames(allvars)==rownames(matches)) < nrow(allvars)){
    stop("Verify inputs: rownames of 'allvars' and 'matches' do not match.")
  }

  # Modify matches data frame if secondarymatch = TRUE
  if (secondarymatch){
    matches <- data.frame(x = matches$x, y = matches$y, target_cell = matches$target_cell,
                          subset_cell = matches$subset_cell_secondary,
                          matching_quality = matches[,quality_name])
  }

  # Create new allvars dataframe with matches column
  allvars1 <- cbind(allvars,subset_cell = as.numeric(matches$subset_cell))

  # Exclude poor matches
  matchedonly <- allvars1[matches$matching_quality <= matching_distance, ]
  # Create results dataframe to store standard deviations:
  results <- data.frame(matrix(rep(NA,(ncol(allvars1)-4)*2),ncol = (ncol(allvars1)-4)))
  colnames(results)[1:(ncol(allvars1)-4)] <- colnames(allvars1)[4:(ncol(allvars1)-1)]
  rownames(results) <- c('allcells','matchedcells')

  # Calculate SD of diffs for all and/or for matched cells:
  if (subset_in_target){
  for (i in 4:(ncol(allvars1)-1)){
    results[1,i-3] <- sd(allvars1[as.character(allvars1[,ncol(allvars1)]),i]-allvars1[,i], na.rm = T)
    results[2,i-3] <- sd(matchedonly[as.character(matchedonly[,ncol(matchedonly)]),i]-matchedonly[,i], na.rm = T)
    }
  }else if (!subset_in_target){
    rownames(subsetcells) <- subsetcells[,subsetcells_id]
    for (i in 4:(ncol(allvars1)-1)){
    results[1,i-3] <- sd(subsetcells[as.character(allvars1[,ncol(allvars1)]),colnames(allvars1)[i]]-allvars1[,i], na.rm = T)
    results[2,i-3] <- sd(subsetcells[as.character(matchedonly[,ncol(matchedonly)]),colnames(allvars1)[i]]-matchedonly[,i], na.rm = T)
    }
  }
  # Remove any columns with all NAs
  rmvec <- NA
  for (i in 1:ncol(results)){
    if (sum(is.na(results[,i])) == nrow(results)){
      rmvec <- append(rmvec,i)
    }
  }
  if (!is.na(rmvec[1])){
  results <- results[,-rmvec]
  }
  # Plot if desired
  if (plot_diffs){
  sd_barplot(results)
  }
  return(results)
  }



#' Internal function for \code{\link{evaluateMatching}}
#'
#' Plots a horizontal barplot of output from \code{\link{evaluateMatching}}
#'
#'
#' @param results data frame. Output from \code{\link{evaluateMatching}}
#'
#'
#' @return a horizonal barplot of the standard deviation of differences between
#' target and matched subset cells.
#'
#'
#' @examples
#'sd_barplot(results)

sd_barplot <- function(results){
  # Calculate xmax_val
  xmax_val <- max(results, na.rm = T)*1.1
  # Create barplot showing standard deviation of differences between Target and Subset cells
  par(mar = c(2,max(nchar(names(results)))/3,2,1), mgp = c(1.5,0.2,0), tcl = 0.3,
      lwd =1, mfrow = c(1,1))
  barplot(as.matrix(results[,ncol(results):1]), beside = T, horiz = T, col = rep(c("grey",0),19),
          names.arg = rev(colnames(results)), las = 1, cex.names = 0.7,
          main = "Standard deviation of differences", xlim = c(0,xmax_val*1.1))
  legend("bottomright", legend = c("Matched only","All cells"), fill= c(0,"grey"), bty = "n")
  box()
}


#' Evaluate matching with geographic distances
#'
#' Calculate 1) distance between target and matched subset cells and 2) distance
#' between the Subset cells matched to each Target cell and the Subset cells
#' matched to the eight adjacent neighbors of that Target cell.
#'
#' @param matches data frame output from the \code{\link{multivarmatch}}
#' function.
#'
#' @param subsetcells if `subset_in_target` is TRUE, this should be a data frame
#' of coordinates (expects coordinates in columns named 'x' and 'y') for Subset
#' cells. May be extracted from output from \code{\link{kpoints}} function or
#' provided separately. Row names should be unique identifiers for each point
#' (unique means no repeats in rownames of subsetcells if `subset_in_target` is
#' TRUE). If `subset_in_target` is FALSE, this should be a data frame of subset
#' cells with column names corresponding exactly to those in `matchingvars` and
#' row names should be unique identifiers (unique means no repeats among all
#' row names in targetcells and matchingvars if `subset_in_target` is FALSE).
#' See `subset_in_target`.
#'
#' @param matchingvars_id character or numeric. Refers to the column in
#' `matchingvars`that provides the unique identifiers for target cells. Defaults
#' to "cellnumbers", which is the unique ID column created by \code{\link{makeInputdata}}.
#'
#' @param subsetcells_id character or numeric, but must be composed of numbers
#' and convertable to numeric. Refers to the column in `subsetcells`that provides
#' the unique identifiers for Subset cells. When `subset_in_target` is TRUE,
#' these ids must be unique from `matchingvars_ids`. Note that if there are
#' repeats between the`matchingvars_id`s and the `subsetcells_id`s, you can paste
#' "00" before the `subsetcells_id`s to ensure they are unique from the
#' `matchingvars_id`s. Defaults to NULL.
#'
#' @param matches data frame output from the \code{\link{multivarmatch}}
#' function.
#'
#' @param quality_name character. Name of the column in the `matches` data frame
#' that contains the matching quality variable to use to evaluate matching
#' "matching_quality" or "matching_quality_secondary".
#' Defaults to "matching_quality"
#'
#' @param subset_in_target boolean. Indicates if Subset cells have been selected
#' from Target cells using \code{\link{kpoints}} function
#'
#' @param matching_distance numeric. Gives the maximum allowable matching quality
#' value (weighted Euclidean distance) between Target and Subset cells. Default
#' value is 1.5.
#'
#' @param longlat boolean. Pass to function in \code{\link[raster-pointDistance]{raster::pointDistance()}}.
#' Indicates if the coordinates are in longitude and latitude format for calculating
#' distances between points. Default value is TRUE and coordinates need to be
#' provided in this format.
#'
#' @param map_distances boolean. Indicates whether a map of distances between
#' Target and matched Subset cells should be plotted. Defaults to TRUE.
#'
#' @param map_neighbor_distances boolean. Indicates whether a map of average
#' distance between the Subset cells matched to each Target cell and the Subset
#' cells matched to the eight adjacent neighbors of that Target cell. Defaults to
#' TRUE.
#'
#' @param which_distance character. One of 'both', 'simple', or 'neighbor'.
#' Determines which distance(s) will be calculated. 'simple' will calculate the
#' dstance between target and matched subset cells, 'neighbor' will calculate the
#' distance between the Subset cells matched to each Target cell and the Subset
#' cells matched to the eight adjacent neighbors of that Target cell. 'both' will
#' calculate both simple and neighbor distances.
#'
#' @param saverasters boolean. Indicates whether to save rasters of the calculated
#' distance metrics. Defaults to FALSE.
#'
#' @param filepath provides path for location where raster will be saved. Defaults
#' to working directory.
#'
#' @param raster_template one of the raster layers used for input data.
#'
#' @return Data frame with the distance between Target and matched Subset cells
#' ('target_to_subset_distance') and the average distance between the Subset cell
#' matched to each Target cell and the Subset cells matched to the eight adjacent
#' Target cells ('avgdistance_to_neighbors'). The first column and the rownames
#' correspond to the unique identifiers for the Target cells, and columns 2 and 3
#' correspond to the 'x' and 'y' coordinates of the Target cells.
#'
#' @examples
#' # Load targetcells data for Target Cells
#' data(targetcells)
#'
#' # Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#'
#' #' # Create vector of matching criteria
#' criteria <- c(0.7,42,3.3,66,5.4,18.4)
#'
#'  # Find solution for k = 200
#' # Note: n_starts should be >= 10, it is 1 here to reduce run time.
#' results1 <- kpoints(matchingvars,criteria = criteria,klist = 200,
#' n_starts = 1,min_area = 50,iter = 50,raster_template = targetcells[[1]])
#'
#'
#' ###################################
#' # First an example where subset_in_target = TRUE
#' # Get points from solution to kpoints algorithm
#' subsetcells <- results1$solutions[[1]]
#'
#' # Find matches and calculate matching quality
#' quals <- multivarmatch(matchingvars, subsetcells, matchingvars_id = "cellnumbers", raster_template = targetcells[[1]], subset_in_target = TRUE)
#'
#' # Look at geographic distances
#' geodist <- evaluateGeoDist(matches = quals, subsetcells = subsetcells,
#'                            subset_in_target = TRUE,
#'                            quality_name = "matching_quality",
#'                            exclude_poor_matches = TRUE,
#'                            matching_distance = 1.5,
#'                            longlat = TRUE, raster_template = targetcells[[1]])
#'
#'
#' ###################################
#' # Now an example where subset_in_target is FALSE
#' # Get points from solution to kpoints algorithm
#' data(subsetcells)
#'
#' # Remove duplicates (representing cells with same climate but different
#' soils--we want to match on climate only)
#' subsetcells <- subsetcells[!duplicated(subsetcells$site_id),]
#'
#' # Pull out matching variables only, with site_id that identifies unique climate
#' subsetcells <- subsetcells[,c("site_id","X_WGS84","Y_WGS84","bioclim_01",
#' "bioclim_04","bioclim_09","bioclim_12",
#' "bioclim_15","bioclim_18")]
#'
#' # Ensure that site_id will be values unique to subsetcells
#' subsetcells$site_id <- paste0("00",subsetcells$site_id)
#'
#' # Find matches and calculate matching quality
#' quals <- multivarmatch(matchingvars, subsetcells=subsetcells,
#' matchingvars_id = "cellnumbers",subsetcells_id = "site_id",
#'                          raster_template = targetcells[[1]],
#'                          subset_in_target = FALSE)
#'
#' # Look at geographic distances
#' geodist <- evaluateGeoDist(matches = quals, subsetcells = subsetcells,
#'                            subsetcells_id = 'site_id',
#'                            subset_in_target = FALSE,
#'                            exclude_poor_matches = TRUE,
#'                            matching_distance = 1.5,
#'                            longlat = TRUE, quality_name = "matching_quality",
#'                            raster_template = targetcells[[1]])

evaluateGeoDist <- function(matches, subsetcells, subsetcells_id = 'site_id',
                            subset_in_target = TRUE,quality_name = "matching_quality",
                            exclude_poor_matches = TRUE,
                            matching_distance = 1.5, longlat = T, raster_template = NULL,
                            map_distances = TRUE, map_neighbor_distances = TRUE,
                            which_distance = "both",saverasters = FALSE,
                            filepath = getwd(), overwrite = FALSE){
  if (which_distance == "simple" | which_distance == "both"){
  # Create matrix of lat/long for cells with matching quality <= 1.5
  # First two columns are coordinates of matched Subset cells
  # Second two columns are coordinates of Target cells
  if (subset_in_target){
    if (exclude_poor_matches){
       pts1 <- cbind(matches[as.character(matches$subset_cell),][matches$matching_quality <= matching_distance,c(1,2)],
                    matches[matches$matching_quality <= matching_distance,c(1,2)])
       rownames(pts1) <- matches[matches$matching_quality <= matching_distance,]$target_cell
    } else if (!exclude_poor_matches){
      pts1 <- cbind(matches[as.character(matches$subset_cell),c(1,2)],
                    matches[,c(1,2)])
      rownames(pts1) <- matches$target_cell
    }
  } else if (!subset_in_target){
    rownames(subsetcells) = subsetcells[,subsetcells_id]
    if (exclude_poor_matches){
      pts1 <- cbind(subsetcells[as.character(as.numeric(matches$subset_cell)),][matches$matching_quality <= matching_distance,c(2,3)],
                    matches[matches$matching_quality <= matching_distance,c(1,2)])
      rownames(pts1) <- matches[matches$matching_quality <= matching_distance,]$target_cell
    } else if (!exclude_poor_matches){
      pts1 <- cbind(subsetcells[as.character(as.numeric(matches$subset_cell)),c(2,3)],
                    matches[,c(1,2)])
      rownames(pts1) <- matches$target_cell
    }
  }


  ## Step 1: Calculate distances between Target and matched Subset cells
  # Calculate distances between pairs of points
  alldistsqual <- raster::pointDistance(p1 = pts1[,1:2], p2 = pts1[,3:4],lonlat = longlat)

  # Convert to km
  distkm <- alldistsqual/1000

  #Create map of distances:
  distmatrix <- cbind(pts1[,3:4], distkm)
  colnames(distmatrix)[3] <- "distance"

  # Create spatial points dataframe from differences:
  ptsx <- sp::SpatialPointsDataFrame(distmatrix[,1:2],
                                     data = data.frame(dist = distkm),
                                     proj4string = raster::crs(raster_template))

  # Rasterize using wydry as a template
  r <- raster::rasterize(ptsx, raster_template, field = ptsx$dist, fun = mean)

  # Save raster, if desired
  if (saverasters){
    raster::writeRaster(r,paste0(filepath,"/Distance_km_subset_to_target.tif"),
                        overwrite = overwrite)
  }

  # Make map of raster, if desired
  if (map_distances){
  # Designate colors
  cols <- c("#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#993404","#662506")

  # set breaks
  bks <- c(0, 6, 12, 18, 36, 72, 144, 288, round(max(distkm)))

  # Plot map:
  legendPlot(r, bks = bks, cols = cols, thisVariable = "Distance to Subset cell (km)")
  }

  # Make results df
  results_simple <- data.frame(target_cell = rownames(distmatrix), x = distmatrix[,1],
                        y = distmatrix[,2], target_to_subset_distance = distmatrix[,3])

  }
  if (which_distance == "neighbors" | which_distance == "both"){
  ## Step 2: Calculate distance between Subset cell matched to each Target cell
    # and the Subset cells matched to that Target cell's 8 adjacent neighbors

  # Make points matrix where cols 3&4 are coordinates of subset cell and
    # 1&2 of target cell
  if (subset_in_target){
    if (exclude_poor_matches){
      pts3 <- cbind(matches[matches[,quality_name] <= matching_distance,c(1,2)],
                    matches[as.character(as.numeric(matches$subset_cell)),][matches[,quality_name] <= matching_distance,c(1,2)])
    } else if (!exclude_poor_matches){
      pts3 <- cbind(matches[,c(1,2)],
                    matches[as.character(as.numeric(matches$subset_cell)),c(1,2)])
    }
  } else if (!subset_in_target){
    rownames(subsetcells) = subsetcells[,subsetcells_id]
    if (exclude_poor_matches){
      pts3 <- cbind(matches[matches[,quality_name] <= matching_distance,c(1,2)],
                    subsetcells[as.character(as.numeric(matches$subset_cell)),][matches[,quality_name] <= matching_distance,c(2,3)])
    } else if (!exclude_poor_matches){
      pts3 <- cbind(matches[,c(1,2)],
                    subsetcells[as.character(as.numeric(matches$subset_cell)),c(2,3)])
    }
  }

  # Create points dataframe
  ptsx <- sp::SpatialPointsDataFrame(pts3[,1:2],
                                     data = data.frame(thispt = rep(1,nrow(pts3))),
                                     proj4string = raster::crs(raster_template))

  # Rasterize using raster_template
  r <- raster::rasterize(ptsx, raster_template, field = 1 , fun = 'first')

  # Use adjacent function to get up to 8 adjacent cells
  # First column gives the focal cell, second column gives any non-NA cell it is adjacent to.
  a <- raster::adjacent(r, cells=as.numeric(row.names(ptsx)), directions=8, pairs=TRUE)

  # Remove nonexistant neighbors (i.e., neighbors not included in marea):
  a1 <- a[which(as.character(a[,2]) %in% row.names(ptsx)),]

  # Calculate distances between subset cells matched to each set of neighbors
  dista <- raster::pointDistance(p1 = pts3[as.character(a1[,1]),3:4], p2 = pts3[as.character(a1[,2]),3:4],lonlat = longlat)

  # Calculate average distance between subset cell matched to target cell and subset cells matched to its 8 neighbors:
  avgdist <- tapply(dista, as.factor(a1[,1]), FUN = mean)

  # Look for cells without neighbors (just FYI):
  nrow(pts3) - length(avgdist)

  # Create distances dataframe:
  adjdist <- data.frame(cellnumber = names(avgdist), avgdist_km = avgdist/1000)
  # Add in coordinates
  adjdist$x <- pts3[adjdist$cellnumber, 1]
  adjdist$y <- pts3[adjdist$cellnumber, 2]

  # Create spatial points dataframe from differences:
  ptsx <- sp::SpatialPointsDataFrame(adjdist[,3:4],
                                     data = data.frame(dist = adjdist[,2]),
                                     proj4string = raster::crs(raster_template))

  # Rasterize using wydry as template
  r <- raster::rasterize(ptsx, raster_template, field = ptsx$dist, fun = mean)

  # Save raster, if desired
  if (saverasters){
    raster::writeRaster(r,paste0(filepath,"/Distance_km_avg_matchedsubset_to_matchedneighbors.tif"),
                        overwrite = overwrite)
  }

  # Make map of raster, if desired
  if (map_distances){
    # Designate colors
    cols <- c("#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#993404","#662506")

    # set breaks
    bks <- c(0, 6, 12, 18, 36, 72, 144, 288, ceiling(max(distkm, na.rm = T)))

    # Make figure
    legendPlot(r, bks = bks, cols = cols, thisVariable = "Mean neighbor distance (km)")
  }
  # Create results df
  results_neighbors <- adjdist
  colnames(results_neighbors)[1] <- "target_cell"
  }
  if (which_distance == "both"){
    results <- merge(results_simple, adjdist[,c(1,2)], by.x = "target_cell",
                     by.y = "cellnumber", all.x = T, all.y = T)
    results <- results[order(as.numeric(results$target_cell)),]
    rownames(results) <- results$target_cell
    colnames(results)[5] <- "avgdistance_to_neighbors"
  } else if (which_distance == "simple"){
    results <- results_simple
    rownames(results) <- results$target_cell
  } else if (which_distance == "neighbors"){
    results <- adjdist[,c(1,3,4,2)]
    colnames(results)[c(1,4)] <- c("target_cell", "avgdistance_to_neighbors")
  }

  return(results)
  }
