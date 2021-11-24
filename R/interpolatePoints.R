#' Interpolate output results
#'
#' Use multivariate matching to interpolate simulation output results to high
#' resolution maps
#'
#'
#' @param matches data frame. Output from `multivarmatch` or `secondaryMatching`
#' functions.
#'
#' @param output_results data frame. Simulation output results for all simulated
#' sites (Subset cells). The first column and the rownames should correspond to
#' the unique identifiers for the Subsetcells. Importantly, these identifiers
#' need to match the identifiers in the 'subset_cell' column of the 'matches'
#' data frame.
#'
#' @param exclude_poor_matches boolean. Indicates whether poor matches (defined
#' as Target cells that are more than the designated 'matching_distance' from
#' their matched Subset cell) should be excluded from matching. Defaults to TRUE.
#'
#' @param subset_cell_names character. This is the name of the column in the
#' 'matches' data frame that provides the unique identity of the Subset cells
#' matched to each Target cell. Defaults to "subset_cell". When 'matches' is the
#' output from `secondaryMatching`, this should be 'subset_cell_secondary'.
#'
#' @param matching_quality_name character. This is the name of the column in the
#' 'matches' data frame that provides the matching quality between the Subset cells
#' and Target cells. Defaults to ""matching_quality"". When 'matches' is the
#' output from `secondaryMatching`, this should be 'matching_quality_secondary'.
#'
#' @param matching_distance numeric. Gives the maximum allowable matching quality
#' value (weighted Euclidean distance) between Target and Subset cells. Default
#' value is 1.5.
#'
#' @param raster_template one of the raster layers used for input data.
#'
#' @param plotraster boolean. Indicates if raster should be plotted to a map.
#' Defaults to TRUE.
#'
#' @param filepath provides path for location where raster will be saved. Defaults
#' to working directory.
#'
#' @param overwrite boolean. Indicates whether saving the rasters will be allowed
#' to overwrite existing files with the same name. Defaults to FALSE.
#'
#' @return raster files of interpolated output variables.
#'
#'
#' @examples
#' # Load targetcells data for Target Cells (from rMultivariateMatchingAlgorithms package)
#' data(targetcells)
#' # Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#' # Subset to include only matching variables
#' matchingvars <- allvars[,c("cellnumbers","x","y","bioclim_01","bioclim_04",
#'  "bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
#' # Create vector of matching criteria
#' criteria <- c(0.7,42,3.3,66,5.4,18.4)
#'
#'
#' # For an example where subset cells were generated from kpoints function
#' # Find solution for k = 200
#' results1 <- kpoints(matchingvars,criteria = criteria,klist = 200,n_starts = 1,
#'                     min_area = 50,iter = 50,raster_template = targetcells[[1]])
#' # Get points from solution to kpoints algorithm
#' subsetcells <- results1$solutions[[1]]
#' # Find matches and calculate matching quality
#' quals <- multivarmatch(matchingvars, subsetcells, criteria = criteria,
#'                          matchingvars_id = "cellnumbers",
#'                          raster_template = targetcells[[1]],
#'                          subset_in_target = TRUE)
#' # Create toy data set of "output variables"
#' # There are really just climate variables from the 'targetcells' rasters,
#' # but we will treat them as output variables to illustrate the method
#' output_results <- allvars[rownames(subsetcells),c("cellnumbers","bioclim_02","bioclim_03","bioclim_16","bioclim_17")]

#' # Interpolate simulation output to rasters
#' interpolatePoints(matches = quals, output_results = output_results,
#'                   exclude_poor_matches = TRUE,
#'                   subset_cell_names = "subset_cell",
#'                   matching_quality_name = "matching_quality",
#'                   matching_distance = 1.5, raster_template = targetcells[[1]],
#'                   plotraster = TRUE, filepath = getwd(),
#'                   overwrite = FALSE)
#'
#'
#'
# # For example where subset cells were not generated from kpoint function
# # Get points from solution to kpoints algorithm
# data(subsetcells)
#' # Pull results from subsetcells
#' output_results <- subsetcells[,c("site_ids","Dryprop","CwetWinter","CdrySummer",
#'                                  "Cwet8","Dryall","Dryany")]
#' rownames(output_results) <- paste0("00",output_results$site_ids)
#'
#' # Remove duplicates (representing cells with same climate but different soils--
#' # we want to match on climate only)
#' subsetcells <- subsetcells[!duplicated(subsetcells$site_id),]
#' # Pull out matching variables only, with site_id that identifies unique climate
#' subsetcells1 <- subsetcells[,c("site_id","X_WGS84","Y_WGS84",names(matchingvars)[4:9])]
#' # Ensure that site_id will be values unique to subsetcells
#' subsetcells1$site_id <- paste0("00",subsetcells$site_id)
#' # Find matches and calculate matching quality
#' quals <- multivarmatch(matchingvars, subsetcells=subsetcells1,
#'                          criteria = criteria,
#'                          matchingvars_id = "cellnumbers",
#'                          subsetcells_id = "site_id",
#'                          raster_templat = targetcells[[1]],
#'                          subset_in_target = FALSE)
#'
#' # Bring in matching variables for secondary matching
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
#'
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
#'
#' # Interpolate simulation output to rasters
#' interpolatePoints(matches = quals2, output_results = output_results,
#'                   exclude_poor_matches = FALSE,
#'                   subset_cell_names = "subset_cell_secondary",
#'                   matching_quality_name = "matching_quality_secondary",
#'                   matching_distance = 1.5, raster_template = targetcells[[1]],
#'                   plotraster = TRUE, filepath = getwd(),
#'                   overwrite = FALSE)


interpolatePoints <- function(matches = NULL, output_results = NULL,
                              exclude_poor_matches = TRUE,
                              subset_cell_names = "subset_cell",
                              matching_quality_name = "matching_quality",
                              matching_distance = 1.5,
                              raster_template = NULL,
                              plotraster = TRUE,
                              filepath = getwd(),
                              overwrite = FALSE){
  if (exclude_poor_matches){
    matches <- matches[matches[,matching_quality_name] <= matching_distance,]
  }

  # Interpolate variables
  for (i in 2:ncol(output_results)){
    print(paste0("Now interpolating ", names(output_results)[i],"."))
  theseresults <- output_results[matches[,subset_cell_names],i]
  # Create spatial points dataframe from differences:
  ptsx <- sp::SpatialPointsDataFrame(matches[,1:2],
                                     data = data.frame(results = theseresults),
                                     proj4string = raster::crs(raster_template))
  # Rasterize
  r <- raster::rasterize(ptsx, raster_template, field = "results",
                         fun = mean, na.rm = T)
  names(r) <- names(output_results[i])
  if (plotraster){
      if (nchar(as.character(floor(mean(theseresults)))) > 1){
        rounding = 0
      } else {
        rounding = 2
      }
    legendPlot(r, thisVariable = names(output_results)[i], round_dec = rounding)
  }
  raster::writeRaster(r,paste0(filepath,"/",names(output_results)[i],".tif"),
                      overwrite = overwrite)
  }
}
