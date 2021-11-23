#' Run a second level of matching
#'
#' In a case where matching on different kinds of variables is required, it can
#' be useful to run a second level of matching. This function is only useful for
#' a case where experimental (simulation) design includes multiple "treatments"
#' (e.g., soil types, aspect) for each set of environmental characteristics
#' (e.g., climate).
#'
#'
#' @param secondaryvars data frame generated using \code{\link{makeInputdata}},
#' or a subset of such a data frame, and/or formatted such that: column 1 and
#' rownames are 'cellnumbers' extracted using the \code{\link{raster::extract}}
#' function, columns 2 and 3 correspond to x and y coordinates, and additional
#' columns correspond to a secondary set of matching variables extracted using
#' the \code{\link{raster::rasterToPoints}} function. These data represent
#' Target cells.
#'
#' @param matches data frame. Output returned from \code{\link{matchingquality}}
#' using primary matching variables (e.g., climate variables).
#'
#' @param subsetcells data frame with columns that correspond to those in
#' `secondaryvars`. Currently, there is no functionality for `subset_in_target` =
#' TRUE, so subset cells should represent a separate set of simulated (Subset)
#' cells. This function is designed to handle experimental design in which there
#' is a "reference" treatment that is site-specific (e.g., site-specific soils)
#' and a series of "other" treatments that are the same across all Subset cells.
#' Thus, the subsetcells data frame should contain secondary variables that
#' correspond to the site-specific treatment only.
#'
#' @param secondaryvars_id character or numeric. Refers to the column in
#' `matchingvars`that provides the unique identifiers for Target cells. Defaults
#' to "cellnumbers", which is the unique ID column created by \code{\link{makeInputdata}}.
#'
#' @param reference_treatment character. Designates the reference
#' treatment identifier that will be pasted onto the `subsetcells_id` to generate
#' unique identifiers for each treatment. Default value is "1".
#'
#' @param subsetcells_id character or numeric. Refers to the column in
#' `subsetcells`that provides the unique identifiers for Target cells. Defaults
#' to NULL.
#'
#' @param other_treatments data frame. Provides secondary variables for "other"
#' treatments that are common among all Subset cells (e.g., a set of soil types
#' that are simulated for each site). The column names should correspond to the
#' secondary variable names and the rownames will designate the treatment identifiers
#' that will be pasted onto the `subsetcells_id` to generate unique identifiers
#' for each treatment.
#'
#' @param criteria single value or vector of length equal to the number of
#' secondary variables, where values corresponds to the matching criterion for
#' each secondary variable `secondaryvars`. If a single value, this will be used
#' as matching criteria for all variables. Default value is 1, corresponding to
#' using raw data for matching.
#'
#' @param raster_template one of the raster layers used for input data.
#'
#' @param subset_in_target boolean. Defaults to FALSE. There is no functionality
#' for TRUE at this time.
#'
#' @param is_loocv boolean. Indicates whether the function is being used as part
#' of leave-one-out cross-validation. Usually only passed to this function from
#' \code{\link{loocv}}.
#'
#' @param saveraster boolean. Indicates if raster of matching quality should be
#' saved to file. Default is FALSE
#'
#' @param plotraster boolean. Indicates if raster of matching quality should be
#' saved to file. Default is TRUE.
#'
#' @param filepath provides path for location where raster will be saved. Defaults
#' to working directory.
#'
#' @param ... additional parameters to pass to legendPlot.
#'
#' @return Data frame of Target cells with coordinates ('x','y'), cellnumber of
#' Target cell ('target_cell'), unique id of matched Subset cell ('subset_cell'),
#' matching quality ('matching_quality'), unique id of Subset cell matched with
#' secondary matching criteria ('subset_cell_secondary'), and matching quality of
#' this secondary match ('matching_quality_secondary). Will save a raster of
#' matching quality if `saveraster` is TRUE and plot a map of matching quality if
#' `plotraster` is TRUE.
#'
#' @examples
#' # Load targetcells data for Target Cells (from rMultivariateMatchingAlgorithms package)
#' data(targetcells)
#' #Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#' # Restrict data to matching variables of interest
#' matchingvars <- allvars[,c("cellnumbers","x","y","bioclim_01","bioclim_04","bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
#' # Create vector of matching criteria
#'
#' # For example with subset_in_target = FALSE (subset_in_target = TRUE is not
#' # functional at this time)
#'
#' # Get points from solution to kpoints algorithm
#' data(subsetcells)
#' # Remove duplicates (representing cells with same climate but different soils--
#' # we want to match on climate only)
#' subsetcells <- subsetcells[!duplicated(subsetcells$site_id),]
#' # Pull out matching variables only, with site_id that identifies unique climate
#' subsetcells <- subsetcells[,c("site_id","X_WGS84","Y_WGS84",names(matchingvars)[4:9])]
#' # Ensure that site_id will be values unique to subsetcells
#' subsetcells$site_id <- paste0("00",subsetcells$site_id)
#' # Find matches and calculate matching quality
#' quals <- matchingquality(matchingvars, subsetcells=subsetcells, matchingvars_id = "cellnumbers",subsetcells_id = "site_id",
#'                          raster_template = targetcells[[1]], subset_in_target = FALSE)
#'
#' # Subset to include only secondaryvars
#' secondaryvars <- allvars[,c("cellnumbers","x","y","sand","clay")]
#'
#' # Bring in secondary id variable from subsetcells
#' data(subsetcells)
#' # Remove duplicates (keeping only site-specific soils with site_ids ending
#' # in ".1").
#' subsetcells <- subsetcells[!duplicated(subsetcells$site_id),]
#' # Pull out matching variables only, with site_id that identifies unique climate
#' subsetcells <- subsetcells[,c("site_id","X_WGS84","Y_WGS84",
#'                               "sand","clay"),]
#' # Convert sand and clay to percentage from fraction
#' subsetcells$sand <- subsetcells$sand*100
#' subsetcells$clay <- subsetcells$clay*100
#'
#' # Make sure subsetcell ids are unique
#' subsetcells$site_id <- paste0("00",subsetcells$site_id)

#' # Bring in "other" treatments
#' data(setsoiltypes)
#' other_treatments = setsoiltypes
#' # Calculate criteria
#' criteria = c((max(secondaryvars$sand,na.rm = T)-min(secondaryvars$sand,na.rm = T))/10,
#'              (max(secondaryvars$clay,na.rm = T)-min(secondaryvars$clay,na.rm = T))/10)
#'
#' # Run secondary matching on soils data
#' quals2 <- secondaryMatching(secondaryvars = secondaryvars, matches = quals,
#'                             subsetcells=subsetcells,subsetcells_id = "site_id",
#'                             subset_in_target = FALSE, criteria = criteria,
#'                             raster_template = targetcells[[1]],
#'                             reference_treatment = "1",
#'                             other_treatments = other_treatments)



secondaryMatching <- function(secondaryvars = NULL, matches = NULL, subsetcells = NULL,
                              secondaryvars_id = "cellnumbers", reference_treatment = "1",
                              subsetcells_id = NULL, criteria = 1,
                              other_treatments = NULL, is_loocv = FALSE,
                              raster_template,subset_in_target = FALSE,
                              saveraster=FALSE,plotraster=TRUE,
                              filepath=getwd(), ...){
  # Check subset_in_target
  if (subset_in_target){
    stop("'subset_in_target = TRUE' does not have functionality at this time.")
  }

  # If no standardization or if standardization is the same for all matching variables
  if (length(criteria)==1){
    criteria = rep(criteria, (ncol(secondaryvars)-3))
  }
  # If incorrect number of matching criteria provided, stop fuunction.
  if (length(criteria) > 1 && length(criteria) != (ncol(secondaryvars)-3)){
    stop("Number of matching criteria unequal to number of matching variables; cannot complete matching.")
  }

  # verify secondaryvars data frame has proper setup:
  if (is.null(secondaryvars$x) | is.null(secondaryvars$y)) {
    stop("Coordinates not found as columns 'x' and 'y' in secondaryvars data frame.")
  }
  if (max(which(colnames(secondaryvars) %in% c('x','y'))) > 3){
    stop("Coordinates not found in correct columns of secondaryvars.")
  }
  if(is.null(secondaryvars_id)){
    stop("secondaryvar_id is not defined")
  }else if (which(colnames(secondaryvars) == secondaryvars_id) > 3 ){
    stop("secondaryvar_id not found in correct column of secondaryvars.")
  }
  for (ci in 4:ncol(secondaryvars)){
    if (!is.numeric(secondaryvars[,ci])){
      stop("Matching variable '",colnames(secondaryvars)[ci],"' is not numeric.")
    }
  }

  # If subset_in_target is FALSE (subset cells provided in separate data frame)
  if (!subset_in_target){
    if (is.null(subsetcells_id)){
      stop("Missing 'subsetcells_id', cannot complete matching.")
    }
    # Coerce id columns to character for comparison:
    secondaryvars[,secondaryvars_id] <- as.character(secondaryvars[,secondaryvars_id])
    subsetcells[,subsetcells_id] <- as.character(subsetcells[,subsetcells_id])

    # Check that id variable columns are all unique
    if (sum(secondaryvars[,secondaryvars_id] %in% subsetcells[,subsetcells_id]) > 0 && !is_loocv){
      warning("Non-unique identifiers in subsetcells and secondaryvars, matching may fail.")
    }

    # Set secondary id columns to rownames
    if (!is_loocv){
    rownames(secondaryvars) <- secondaryvars[,secondaryvars_id]
    rownames(subsetcells) <- subsetcells[,subsetcells_id]
    }
  }

  # Calculate distance between reference treatment and target cells
  site_dist <- (secondaryvars[,4]/criteria[1]-subsetcells[matches$subset_cell,4]/criteria[1])^2
  for (i in 5:ncol(secondaryvars)){
    site_dist = site_dist +
      (secondaryvars[,i]/criteria[i-3]-subsetcells[matches$subset_cell,i]/criteria[i-3])^2
  }
  site_dist = sqrt(site_dist)

  # Set names to the Subset cell name that corresponds to reference treatment
  names(site_dist) <- paste0(matches$subset_cell,".1")

  # Remove missing values
  site_dist <- na.omit(site_dist)

  # Combine reference and other treatment variables
  alltmts <- rbind(other_treatments,secondaryvars[,c(4:ncol(secondaryvars))])

  # Transform alltmts:
  alltmts1 <- alltmts[,1]/criteria[1]
  for (i in 2:ncol(alltmts)){
    alltmts1 <- cbind(alltmts1, alltmts[,i]/criteria[i])
  }
  rownames(alltmts1) <- rownames(alltmts)
  # Remove any missing values
  alltmts1 <- na.omit(alltmts1)

  # Calculate distances between set soil types and matched cells:
  xdist1 <- distances::distances(alltmts1, id_variable = row.names(alltmts1))
  cell_numbers <- rownames(alltmts1)
  # Find nearest neighbor from among other treatments
  neighbors_tmts <-
    distances::nearest_neighbor_search(xdist1, k = 1,
                            search_indices = c(1:nrow(other_treatments)),
                            query_indices =
                              c((nrow(other_treatments)+1):nrow(alltmts1)))

  # Create dataframe with matched site variables corresponding to each cell
  site_vars2a <- alltmts1[neighbors_tmts,]

  # Calculate distance between other treatments and matched cells
  other_dist <- (site_vars2a[,1]-na.omit(secondaryvars[,4])/criteria[1])^2
  for (i in 5:ncol(secondaryvars)){
    other_dist <- other_dist +(site_vars2a[,i-3]-na.omit(secondaryvars[,i])/criteria[5-3])^2
  }
  names(other_dist) <- paste0(matches$subset_cell[complete.cases(secondaryvars)],".",names(other_dist))

  # Determine min.dist between reference treatment distance and nearest other
  # treatments distance
  min_dist <- apply(cbind(site_dist, other_dist), 1, FUN = which.min)

  # Create names matrix to determine best match
  names_matrix <- cbind(names(site_dist), names(other_dist))
  names_selection_matrix <- cbind(1:length(min_dist), min_dist)
  secondary_matches <- names_matrix[names_selection_matrix]

  # Calculate secondary matching quality
  secondary_matching <- cbind(site_dist,other_dist)[names_selection_matrix]

  # Remove any site without soils information from matches
  matches <- matches[complete.cases(secondaryvars),]
  # Create a matches column for new matches
  matches$subset_cell_secondary = secondary_matches
  # Create a matching quality column
  matches$matching_quality_secondary = secondary_matching

  # Add in matches for other treatments if is_loocv
  if (is_loocv){
    matches$tmp_target_cell <- paste0(matches$target_cell,".",reference_treatment)
    rownames(matches) <- matches$tmp_target_cell
    # Now loop through and add other treatments:
    matches1 <- matches
    for (i in rownames(other_treatments)){
    thesematches <- matches
    thesematches$subset_cell_secondary <- paste0(thesematches$subset_cell,".",i)
    thesematches$tmp_target_cell <- paste0(thesematches$target_cell,".",i)
    thesematches$matching_quality_secondary <- 0
    row.names(thesematches) <- thesematches$tmp_target_cell
    matches1 <- rbind(matches1, thesematches)
    }
    matches1$target_cell <- matches1$tmp_target_cell
    matches <- matches1[,-which(colnames(matches1) == "tmp_target_cell")]
  }

  if (saveraster | plotraster){
  # Create spatial points dataframe
  ptsx <- sp::SpatialPointsDataFrame(matches[,1:2],
                                     data = matches,
                                     proj4string = crs(raster_template))

  # Create raster of qual (matching quality) using wydry as a template
  r <- raster::rasterize(ptsx, raster_template,
                         field = matches$matching_quality_secondary, fun = mean)

  if (saveraster){
    raster::writeRaster(r,paste0(filepath,"/Matchingquality_secondary.tif"))
  }

  if (plotraster){
    # Designate colors, breaks, and plotting settings
    cols = rev(c("#d7191c","#fdae61","#abd9e9","#2c7bb6"))
    bks = c(0,0.5,1,1.5,5)
    thisVariable = "Matching quality"
    legendPlot(r, bks = bks, cols = cols, thisVariable = "Secondary matching quality",matchingQ = TRUE, ...)
  }
  }
  return(matches)
}
