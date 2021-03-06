% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluateMatching.R
\name{multivarmatch}
\alias{multivarmatch}
\title{Identify matches, calculate and map matching quality}
\usage{
multivarmatch(
  matchingvars = NULL,
  subsetcells = NULL,
  matchingvars_id = "cellnumbers",
  subsetcells_id = NULL,
  criteria = 1,
  n_neighbors = 1,
  raster_template,
  subset_in_target = TRUE,
  saveraster = FALSE,
  plotraster = TRUE,
  addpoints = FALSE,
  filepath = getwd(),
  overwrite = FALSE,
  ...
)
}
\arguments{
\item{matchingvars}{data frame generated using \code{\link{makeInputdata}} or
formatted such that: column 1 and rownames are 'cellnumbers' extracted using the
\code{\link[raster]{extract}} function, columns 2 and 3 correspond to x and y
coordinates, and additional columns correspond to potential matching variables
extracted using the \code{\link[raster]{rasterToPoints}} function. These data
represent Target cells.}

\item{subsetcells}{if \code{subset_in_target} is TRUE, this should be a data frame
of coordinates (expects coordinates in columns named 'x' and 'y') for Subset
cells. May be extracted from output from \code{\link{kpoints}} function or
provided separately. Row names should be unique identifiers for each point
(unique means no repeats in rownames of subsetcells if \code{subset_in_target} is
TRUE). If \code{subset_in_target} is FALSE, this should be a data frame of subset
cells with column names corresponding exactly to those in \code{matchingvars} and
row names should be unique identifiers (unique means no repeats among all
row names in targetcells and matchingvars if \code{subset_in_target} is FALSE).
See \code{subset_in_target}.}

\item{matchingvars_id}{character or numeric. Refers to the column in
\code{matchingvars}that provides the unique identifiers for Target cells. Defaults
to "cellnumbers", which is the unique ID column created by \code{\link{makeInputdata}}.}

\item{subsetcells_id}{character or numeric, but must be composed of numbers
and convertable to numeric. Refers to the column in \code{subsetcells}that provides
the unique identifiers for Subset cells. When \code{subset_in_target} is TRUE,
these ids must be unique from \code{matchingvars_ids}. Note that if there are
repeats between the\code{matchingvars_id}s and the \code{subsetcells_id}s, you can paste
"00" before the \code{subsetcells_id}s to ensure they are unique from the
\code{matchingvars_id}s. Defaults to NULL.}

\item{criteria}{single value or vector of length equal to the number of matching variables,
where values corresponds to the matching criterion for each matching variable
in x. If a single value, this will be used as matching criteria for all variables.
Default value is 1, corresponding to using raw data for matching.}

\item{n_neighbors}{numeric. The number of neighbors to search for in matching.
Default value is 1 and this setting should be used for matching. Option for 2+
is only included for leave-one-out cross-validation \code{\link{loocv}}.}

\item{raster_template}{one of the raster layers used for input data.}

\item{subset_in_target}{boolean. Indicates if Subset cells have been selected
from Target cells using \code{\link{kpoints}} function.}

\item{saveraster}{boolean. Indicates if raster of matching quality should be
saved to file.}

\item{plotraster}{boolean. Indicates if raster should be plotted to a map.}

\item{addpoints}{boolean. Indicates if Subset cells should be added to the plot
as points. Defaults to FALSE.}

\item{filepath}{provides path for location where raster will be saved. Defaults
to working directory.}

\item{overwrite}{boolean. Indicates whether \code{\link[raster]{writeRaster}}
should overwrite existing files with the same name in \code{filepath}. Defaults to
FALSE.}

\item{...}{additional arguments to pass to \code{\link{legendPlot}}.}
}
\value{
Data frame of Target cells with coordinates ('x','y'), cellnumber of
Target cell ('target_cell'), unique id of matched Subset cell ('subset_cell')
and matching quality ('matching_quality'). Will save a raster of matching
quality if \code{saveraster} is TRUE and plot a map of matching quality if
\code{plotraster} is TRUE.
}
\description{
Identifies matches from Subset cells for all Target cells, then calculates
matching quality (weighted Euclidean distance between Target and matched
Subset cells), with options to plot map of matching quality and save a raster
of matching quailty.
}
\examples{
# Load targetcells data for Target Cells
data(targetcells)

# Create data frame of potential matching variables for Target Cells
allvars <- makeInputdata(targetcells)

# Restrict data to matching variables of interest
matchingvars <- allvars[,c("cellnumbers","x","y","bioclim_01","bioclim_04",
                       "bioclim_09","bioclim_12","bioclim_15","bioclim_18")]

# Create raster_template
raster_template <- targetcells[[1]]

# Create vector of matching criteria
criteria <- c(0.7,42,3.3,66,5.4,18.4)

# Find solution for k = 200
# Note: n_starts should be >= 10, it is 1 here to reduce run time.
results1 <- kpoints(matchingvars,criteria = criteria,klist = 200,
                    n_starts = 1,min_area = 50,iter = 50,
                    raster_template = raster_template)

###################################
# First an example where subset_in_target = TRUE

# Get points from solution to kpoints algorithm
subsetcells <- results1$solutions[[1]]

# Find matches and calculate matching quality
quals <- multivarmatch(matchingvars, subsetcells, criteria = criteria,
                        matchingvars_id = "cellnumbers", addpoints = FALSE,
                        raster_template = raster_template,
                        subset_in_target = TRUE)

###################################
# Now an example where subset_in_target is FALSE
# Remove previous subsetcells
rm(subsetcells)

# Get Subset cells data
data(subsetcells)

# Remove duplicates (representing cells with same climate but different soils--
# we want to match on climate only)
subsetcells <- subsetcells[!duplicated(subsetcells$site_id),]

# Pull out matching variables only, with site_id that identifies unique climate
subsetcells <- subsetcells[,c("site_id","X_WGS84","Y_WGS84","bioclim_01",
                           "bioclim_04","bioclim_09","bioclim_12",
                           "bioclim_15","bioclim_18")]

# Ensure that site_id will be values unique to subsetcells
subsetcells$site_id <- paste0("00",subsetcells$site_id)

# Find matches and calculate matching quality
quals <- multivarmatch(matchingvars, subsetcells, criteria = criteria,
                         matchingvars_id = "cellnumbers",
                         subsetcells_id = "site_id",
                         raster_template = raster_template,
                         subset_in_target = FALSE, addpoints = FALSE)
}
\author{
Rachel R. Renne
}
