#' Set up input data for Target cells
#'
#' Prepare data frame of potential matching variables for Target cells from
#' raster files
#'
#' @param x RasterStack that covers the area of interest and includes all potential
#' matching variables. Large datasets may be difficult to handle if memory is
#' limited and it may work better to run the function for several RasterStacks
#' (removing previous RasterStacks before reading in the next one) and combine
#' the output into one final data frame. Note: rasters must all have the same
#' extent, resolution, and crs.
#'
#' @return data frame with a 'cellnumbers' column (these are unique values for each
#' cell in the rasters used as input and will be used as unique identifiers for
#' the Target cells in other functions), x and y coordinates, and values for each
#' of the rasters in the RasterStack. The rownames correspond to the 'cellnumbers'
#' column.
#'
#'
#' @examples
#' # Load targetcells data for Target Cells
#' data(targetcells)
#' # Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#'

# Function to format rasters of potential matching variables for Target Cells
makeInputdata <- function(x){
  if (sum(grep("raster", class(x), ignore.case = TRUE)) < 1){
    stop("Incorrect inputs, x must be a raster.")
  }
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





#' Determine matching criteria
#'
#' Examine summary statistics of matching variables to determine matching criteria
#'
#' @param matchingvars data frame created using \code{\link{makeInputdata}} or
#' formatted such that: rownames are 'cellnumbers' extracted using the
#' \code{\link[raster:extract]{raster::extract()}} function, columns 2 and 3
#' correspond to x and y coordinates, and additional columns correspond to
#' potential matching variables extracted using the
#' \code{\link[raster:rasterToPoints]{raster::rasterToPoints()}} function.
#'
#' @param criteria_list list of matching criteria to test in the \code{\link{kpoints}}
#' function. Each item in the list should be a vector of values (possible matching
#' criteria) corresponding to each matching variable.
#'
#' @param k number of points to use to find solution using \code{\link{kpoints}}
#' function. Default value is 200.
#'
#' @param raster_template one of the raster layers used for input data.
#' See \code{\link[raster:area]{raster::area()}}. Note that 'cellnumbers'
#' column must be present for \code{\link{kpoints}} function to work within this
#' function.
#'
#' @param plot_coverage boolean. Indicates whether the algorithm should display
#' a barplot of the coverage for each set of criteria. Default is TRUE.
#'
#' @param subset_in_target boolean. Indicates if Subset cells have been selected
#' from Target cells using \code{\link{kpoints}} function
#'
#' #' @param subsetcells data frame. Passed to \code{\link{multivarmatch}} function
#' if `subset_in_target` is FALSE. This should be a data frame of subset
#' cells with column names corresponding exactly to those in `matchingvars` and
#' row names should be unique identifiers.
#'
#' @param matchingvars_id character or numeric. Refers to the column in
#' `matchingvars`that provides the unique identifiers for Target cells. Defaults
#' to "cellnumbers", which is the unique ID column created by \code{\link{makeInputdata}}.
#'
#' @param subsetcells_id character or numeric. Refers to the column in
#' `subsetcells`that provides the unique identifiers for Subset cells. Defaults
#' to "site_id".
#'
#' @param matching_distance Gives the maximum allowable matching quality
#' value (weighted Euclidean distance) between Target and Subset cells, when
#' `subset_in_target` is FALSE. Default value is 1 so that output will be
#' comparable to output from `choose_criteria` when `subset_in_target` is TRUE.
#'
#' @param ... accepts additional parameters to \code{\link{kpoints}} function.
#'
#' @return
#'
#' @examples
#' # Load targetcells data for Target Cells
#' data(targetcells)
#' # Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#'
#' # Restrict data to matching variables of interest
#' matchingvars <- allvars[,c("cellnumbers","x","y","bioclim_01","bioclim_04",
#' "bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
#'
#' # Create list of matching criteria to choose:
#' # Look at 2.5%, 5%, & 10% of range and standard deviation for each variable
#' range2.5pct <- apply(matchingvars[,4:ncol(matchingvars)],2,
#' function(x){(max(x)-min(x))*0.025})
#' range5pct <- apply(matchingvars[,4:ncol(matchingvars)],2,
#' function(x){(max(x)-min(x))*0.05})
#' range10pct <- apply(matchingvars[,4:ncol(matchingvars)],2,
#' function(x){(max(x)-min(x))*0.1})
#' stddev <- apply(matchingvars[,4:ncol(matchingvars)],2,sd)
#' criteria_list <- list(range2.5pct, range5pct, range10pct, stddev)
#'
#' ###################################
#' # First an example where subset_in_target = TRUE
#' # Compare coverage with various criteria
#' # Note: n_starts should be >= 10, it is 1 here to reduce run time.
#' results2 <- choose_criteria(matchingvars, criteria_list = criteria_list,
#' n_starts = 1, k = 200, raster_template = targetcells[[1]],
#' subset_in_target = TRUE, plot_coverage = TRUE)
#'
#' ###################################
#' # Now an example where subset_in_target is FALSE
#' # Bring in Subset cell data
#' data(subsetcells)
#' # Remove duplicates (representing cells with same climate but different soils--
#' # we want to match on climate only)
#' subsetcells <- subsetcells[!duplicated(subsetcells$site_id),]
#' # Pull out matching variables only, with site_id that identifies unique climate
#' subsetcells <- subsetcells[,c("site_id","X_WGS84","Y_WGS84"
#'                               names(matchingvars)[4:9])]
#' # Ensure that site_id will be values unique to subsetcells
#' subsetcells$site_id <- paste0("00",subsetcells$site_id)
#' # Run choose_criteria function to evaluate different matching criteria
#' coverage <- choose_criteria(matchingvars = matchingvars,
#'                             criteria_list = criteria_list,
#'                             plot_coverage = TRUE,
#'                             raster_template = targetcells[[1]],
#'                             subset_in_target = FALSE,
#'                             subsetcells = subsetcells,
#'                             matchingvars_id = "cellnumbers",
#'                             subsetcells_id = "site_id")

choose_criteria <- function(matchingvars = NULL, criteria_list = NULL,
                            k = 200,plot_coverage = TRUE,
                            raster_template = NULL,
                            subset_in_target = TRUE,
                            subsetcells = NULL,
                            matchingvars_id = "cellnumbers",
                            subsetcells_id = "site_id",
                            matching_distance = 1,
                            ...){
  if (class(criteria_list) != "list" | is.null(criteria_list)){
    stop("Verify inputs: 'criteria_list' is not a list or is missing.")
  }
  if (subset_in_target){
  criteria_results <- list()
  for (cr in 1:length(criteria_list)){
    print(paste0("Using criteria: ", criteria_list[cr]))
    criteria = criteria_list[[cr]]
    klist = k
    results1 <- kpoints(matchingvars=matchingvars,criteria,klist=klist,...)
    criteria_results$solution_areas[cr] <- results1$solution_areas
    criteria_results$totalarea = results1$totalarea
    criteria_results$k = results1$k
    criteria_results$iter = results1$iter
    criteria_results$n_starts = results1$n_starts
    criteria_results$criteria_list = results1$criteria_list
    criteria_results$min_area = results1$min_area
  }
  if (plot_coverage){
  criteriaplot(criteria_results, criteria_list)
  }
  } else if (!subset_in_target){
    criteria_results <- list()
    areas <- raster::area(raster_template)
    criteria_results[["totalarea"]] <- round(sum(raster::extract(areas, as.numeric(row.names(matchingvars)))))
    for (i in 1:length(criteria_list)){
      # Find matches and calculate matching quality
      quals <- multivarmatch(matchingvars, subsetcells, criteria = criteria_list[[i]],
                             matchingvars_id = matchingvars_id,
                             subsetcells_id = subsetcells_id,
                             raster_template = raster_template,
                             matching_distance = matching_distance,
                             plotraster = FALSE,
                             subset_in_target = FALSE,
                             addpoints = FALSE, ...)

      criteria_results[["solution_areas"]][i] <- round(sum(extract(areas,as.numeric(rownames(quals[quals$matching_quality <= matching_distance,])))))
    }
    if (plot_coverage){
      criteriaplot(criteria_results, criteria_list)
    }
  }
  return(criteria_results)
}


#' Create a barplot of results from \code{\link{choose_criteria}} function
#'
#' Shows a barplot of the proportion of the study area covered by each of the sets
#' of matching criteria that were tested in the \code{\link{choose_criteria}}
#' function.
#'
#'
#' @param x a list of results from the \code{\link{choose_criteria}} function.
#'
#' @param criteria_list. A list of matching criteria tested in the
#' \code{\link{choose_criteria}} function.
#'
#'
#' @return a barplot showing the proportion of the study area represented by
#' different combinations of matching criteria.
#'
#'
#' @examples
#' See \code{\link{choose_criteria}}.


criteriaplot <- function(x, criteria_list){
par(tcl = 0.3, mgp = c(1.5,0.2,0), mar = c(3,3,1,1))
barplot(x[['solution_areas']]/x[['totalarea']],
        main = "Area covered for different criteria",
        xlab = "Criteria (refer to criteria_list)",
        ylab = "Proportion of area covered",
        cex.lab = 1.2, cex.axis = 1.2,
        ylim = c(0,1), names.arg = c(1:length(criteria_list)))
box()
}
