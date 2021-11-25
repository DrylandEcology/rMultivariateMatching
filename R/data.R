#' Climate and soil variables for the state of Wyoming
#'
#' A rasterStack containing 21 rasters of climate and soil variables for drylands
#' in the state of Wyoming at 30-arcsecond resolution. We defined drylands as
#' cells where a 30-arcsecond (~1 km) resolution aridity index (ratio of
#' precipitation to potential evapotranspiration) defined for the 1970-2000 period
#' (Trabucco & Zomer, 2019) was <0.6.
#'
#' The climate variables ("bioclim_XX") correspond to the 19 bioclim variables
#' described by Hijman (2017) in \code{\link[dismo]{biovars}} function (also
#' defined below).
#'
#' @format A rastereStack with 21 layers, each with the following attributes:
#' dimensions : 482, 841, 405362, 21  (nrow, ncol, ncell, nlayers)
#' resolution : 0.008333333, 0.008333333  (x, y)
#' extent     : -111.0583, -104.05, 40.99167, 45.00833  (xmin, xmax, ymin, ymax)
#' crs        : +proj=longlat +datum=WGS84 +no_defs
#'
#' \describe{
#'   \item{bioclim_01}{Mean annual temperature, degrees Celsius}
#'   \item{bioclim_02}{Mean diurnal range (mean of max temp - min temp), degrees Celsius}
#'   \item{bioclim_03}{Isothermality (bio2/bio7) (* 100)}
#'   \item{bioclim_04}{Temperature seasonality (standard deviation *100)}
#'   \item{bioclim_05}{Max temperature of warmest month, degrees Celsius}
#'   \item{bioclim_06}{Min temperature of coldest month, degrees Celsius}
#'   \item{bioclim_07}{Temperature annual range (bio5-bio6), degrees Celsius}
#'   \item{bioclim_08}{Mean temperature of the wettest quarter, degrees Celsius}
#'   \item{bioclim_09}{Mean temperature of driest quarter, degrees Celsius}
#'   \item{bioclim_10}{Mean temperature of warmest quarter, degrees Celsius}
#'   \item{bioclim_11}{Mean temperature of coldest quarter, degrees Celsius}
#'   \item{bioclim_12}{Total (annual) precipitation, mm}
#'   \item{bioclim_13}{Precipitation of wettest month, mm}
#'   \item{bioclim_14}{Precipitation of driest month, mm}
#'   \item{bioclim_15}{Precipitation seasonality (coefficient of variation)}
#'   \item{bioclim_16}{Precipitation of wettest quarter, mm}
#'   \item{bioclim_17}{Precipitation of driest quarter, mm}
#'   \item{bioclim_18}{Precipitation of warmest quarter, mm}
#'   \item{bioclim_19}{Precipitation of coldest quarter, mm}
#'   \item{sand}{Percentage sand content of soil by weight}
#'   \item{clay}{Percentage clay contet of soil by weight}
#' }
#' @source We calculated 30-year (1981 to 2010) normal monthly temperature and
#' precipitation data at 30-arcsecond resolution from DayMet (Thornton et al.,
#' 2018) using Google Earth Engine (Gorelick et al., 2017). We then calculated 19
#' bioclimatic variables using the \code{\link[dismo]{biovars}} function in R package dismo
#' (Hijmans, 2017) We derived soils information for each 30-arcsecond cell from
#' Soilgrids+ products (Hengl et al., 2017). We downloaded 250 m resolution global
#' rasters for sand and clay from 0 to 100 cm, cropped each raster to our study
#' area, aggregated the 250 m cells to ~1 km resolution, and reprojected them to
#' a 30-arcsecond resolution to match the bioclim variables. The Soilgrids+ data
#' provide predicted sand and clay content at seven depths: 0 cm, 5 cm, 15 cm,
#' 30 cm, 60 cm, 100 cm, and 200 cm. We derived site specific soil data by
#' calculating depth-weighted sand and clay content for the top 100 cm using the
#' following formula (from Hengl et al., (2017)):
#'
#' ((5)(L1+L2)+(10)(L2+L3)+(15)(L3+L4)+(30)(L4+L5)+(40)(L5+L6))/200
#'
#' where L1 is the clay or sand content at 0 cm, L2 is the clay or sand content
#' at 5 cm, etc. We then extracted the depth-weighed sand and clay contents for
#' each 30-arcsecond cell. Some cells do not have soils data due to the presence
#' of waterbodies or glaciers.
#'
#' @References
#' Gorelick, N., Hancher, M., Dixon, M.,  Ilyushchenko, S., Thau, D., & Moore, R.
#' (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone.
#' Remote Sens. Environ. 202, pp.18-27. https://doi.org/10.1016/j.rse.2017.06.031
#'
#' Hijmans, R. J., Phillips, S., Leathwick, J., & Elith, J. (2017). dismo: Species
#' Distribution Modeling. R package version 1.1-4. https://CRAN.R-project.org/package=dismo
#'
#' Hengl, T., de Jesus, J. M., Heuvelink, G. B., Gonzalez, M. R., Kilibarda, M.,
#' Blagotić, A., Shangguan, W., Wright, M. N., Geng, X., Bauer-Marschallinger, B.,
#' Guevara, M.A. Vargas, R., MacMillan, R. A., Batjes, N. H., Leenaars, J. G. B.,
#' Ribeiro, E., Wheeler, I., Mantel, S., & Kempen, B. (2017). SoilGrids250m:
#' Global gridded soil information based on machine learning. PLoS one, 12(2),
#' p.e0169748. https://doi.org/10.1371/journal.pone.0169748
#'
#' Thornton, P. E., Thornton, M. M., Mayer, B. W., Wei, Y., Devarakonda, R., Vose,
#' R. S., & Cook, R. B. (2018). Daymet: Daily Surface Weather Data on a 1-km Grid
#' for North America, Version 3. ORNL DAAC, Oak Ridge, Tennessee, USA.
#' https://doi.org/10.3334/ORNLDAAC/1328
#'
#' Trabucco, A., & Zomer, R. (2019). Global Aridity Index and Potential
#' Evapotranspiration (ET0) Climate Database v2. figshare. Fileset.
#' https://doi.org/10.6084/m9.figshare.7504448.v3
"targetcells"



#' Climate, soil, and ecohydrological variables for a 10 x 10 km grid of Wyoming
#'
#' A data frame containing unique ids, geographic coordinates, 19 climate variables,
#' 2 soil variables and 6 ecohydrological variables for drylands in the state of
#' Wyoming at 10-km resolution (Albers Equal Area). We defined drylands as cells
#' where the aridity index (ratio of precipitation to potential evapotranspiration)
#' calculated from SOILWAT2 (Schlaepfer et al., 2012; Schlaepfer & Andrews, 2018;
#' Schlaepfer & Murphy, 2018) was <0.6.
#'
#' The climate variables ("bioclim_XX") correspond to the 19 bioclim variables
#' described by Hijman (2017) in \code{\link[dismo]{biovars}} function (also
#' defined below). Columns 29-34 contain ecohydrological variables derived from
#' SOILWAT2 simulations and described in Bradford et al. (2019). These variables
#' describe the frequency and seasonality of wet (>-1.5MPa) soil conditions within
#' the moisture control section (MCS: soil layers with depth ranging from 10-30 cm
#' for fine textures to 30–90 cm for coarse textures; Soil Survey Staff, 2014).
#'
#' @format A data frame with 10630 rows and 34 variables:
#'
#' \describe{
#'   \item{site_ids}{Unique identifer: site number (5-digits), with soil type
#'   (1-5) designated as a decimal}
#'   \item{Label}{Unique label used to identify sites in SOILWAT2 simulations}
#'   \item{site_id}{Unique identifier: site number (5-digits)}
#'   \item{X_WGS84}{Longitude coordinate in WGS84}
#'   \item{Y_WGS84}{Latitude coordinate in WGS84}
#'   \item{X_AEANAD83}{Longitude coordinate in Albers Equal Area}
#'   \item{Y_AEANAD83}{Latitude coordinate in Alberrs Equal Area}
#'   \item{bioclim_01}{Mean annual temperature, degrees Celsius}
#'   \item{bioclim_02}{Mean diurnal range (mean of max temp - min temp), degrees Celsius}
#'   \item{bioclim_03}{Isothermality (bio2/bio7) (* 100)}
#'   \item{bioclim_04}{Temperature seasonality (standard deviation *100)}
#'   \item{bioclim_05}{Max temperature of warmest month, degrees Celsius}
#'   \item{bioclim_06}{Min temperature of coldest month, degrees Celsius}
#'   \item{bioclim_07}{Temperature annual range (bio5-bio6), degrees Celsius}
#'   \item{bioclim_08}{Mean temperature of the wettest quarter, degrees Celsius}
#'   \item{bioclim_09}{Mean temperature of driest quarter, degrees Celsius}
#'   \item{bioclim_10}{Mean temperature of warmest quarter, degrees Celsius}
#'   \item{bioclim_11}{Mean temperature of coldest quarter, degrees Celsius}
#'   \item{bioclim_12}{Total (annual) precipitation, mm}
#'   \item{bioclim_13}{Precipitation of wettest month, mm}
#'   \item{bioclim_14}{Precipitation of driest month, mm}
#'   \item{bioclim_15}{Precipitation seasonality (coefficient of variation)}
#'   \item{bioclim_16}{Precipitation of wettest quarter, mm}
#'   \item{bioclim_17}{Precipitation of driest quarter, mm}
#'   \item{bioclim_18}{Precipitation of warmest quarter, mm}
#'   \item{bioclim_19}{Precipitation of coldest quarter, mm}
#'   \item{sand}{Percentage sand content of soil by weight}
#'   \item{clay}{Percentage clay contet of soil by weight}
#'   \item{Dryprop}{Proportion of days that all layers within the MCS are dry when
#'    soil temperature at 50 cm >5 degrees Celsius}
#'   \item{CwetWinter}{Number of consecutive days with all MCS layers wet during
#'   the winter}
#'   \item{CdrySummer}{Number of consecutive days with all MCS layers dry during
#'   the summer}
#'   \item{Cwet8}{Number of consecutive days with any layer wet when soil
#'   temperature at 50 cm depth is >8 degrees Celsius}
#'   \item{Dryall}{}
#'   \item{Dryany}{}
#' }
#' @source We calculated 30-year (1981 to 2010) normal monthly temperature and
#' precipitation data at 30-arcsecond resolution from DayMet (Thornton et al.,
#' 2018) using Google Earth Engine (Gorelick et al., 2017). We then calculated 19
#' bioclimatic variables using the \code{\link[dismo]{biovars}} function in R package dismo
#' (Hijmans, 2017) We derived soils information for each 30-arcsecond cell from
#' Soilgrids+ products (Hengl et al., 2017). We downloaded 250 m resolution global
#' rasters for sand and clay from 0 to 100 cm, cropped each raster to our study
#' area, aggregated the 250 m cells to ~1 km resolution, and reprojected them to
#' a 30-arcsecond resolution to match the bioclim variables. The Soilgrids+ data
#' provide predicted sand and clay content at seven depths: 0 cm, 5 cm, 15 cm,
#' 30 cm, 60 cm, 100 cm, and 200 cm. We derived site specific soil data by
#' calculating depth-weighted sand and clay content for the top 100 cm using the
#' following formula (from Hengl et al., (2017)):
#'
#' ((5)(L1+L2)+(10)(L2+L3)+(15)(L3+L4)+(30)(L4+L5)+(40)(L5+L6))/200
#'
#' where L1 is the clay or sand content at 0 cm, L2 is the clay or sand content
#' at 5 cm, etc. We then extracted the depth-weighed sand and clay contents for
#' each 30-arcsecond cell. Some cells do not have soils data due to the presence
#' of waterbodies or glaciers.
#'
#' @Reference
