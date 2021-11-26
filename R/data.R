#' Climate and soil variables for the state of Wyoming
#'
#' A rasterStack containing 21 rasters of climate and soil variables for drylands
#' in the state of Wyoming at 30-arcsecond resolution (WGS84). We defined drylands
#' as cells where a 30-arcsecond (~1 km) resolution aridity index (ratio of
#' precipitation to potential evapotranspiration) defined for the 1970-2000 period
#' (Trabucco & Zomer, 2019) was <0.6.
#'
#' The climate variables ("bioclim_XX") correspond to the 19 bioclim variables
#' described by Hijman (2017) in \code{\link[dismo]{biovars}} function (also
#' defined below).
#'
#' @format A rasterStack with 21 layers, each with the following attributes:\cr
#' \cr
#' **dimensions:** 482, 841, 405362, 21  (nrow, ncol, ncell, nlayers)\cr
#' **resolution:** 0.008333333, 0.008333333  (x, y)\cr
#' **extent:** -111.0583, -104.05, 40.99167, 45.00833  (xmin, xmax, ymin, ymax)\cr
#' **crs:** +proj=longlat +datum=WGS84 +no_defs\cr
#'
#' The 21 layers are described below:
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
#'   \item{clay}{Percentage clay content of soil by weight}
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
#' (5(L1+L2)+10(L2+L3)+15(L3+L4)+30(L4+L5)+40(L5+L6))/200
#'
#' where L1 is the clay or sand content at 0 cm, L2 is the clay or sand content
#' at 5 cm, etc. We then extracted the depth-weighed sand and clay contents for
#' each 30-arcsecond cell. Some cells do not have soils data due to the presence
#' of waterbodies or glaciers.
#'
#' **References**\cr
#' Gorelick, N., Hancher, M., Dixon, M.,  Ilyushchenko, S., Thau, D., & Moore, R.
#' (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone.
#' *Remote Sens. Environ.* 202, pp.18-27. https://doi.org/10.1016/j.rse.2017.06.031
#'
#' Hijmans, R. J., Phillips, S., Leathwick, J., & Elith, J. (2017). dismo: Species
#' Distribution Modeling. R package version 1.1-4. https://CRAN.R-project.org/package=dismo
#'
#' Hengl, T., de Jesus, J. M., Heuvelink, G. B., Gonzalez, M. R., Kilibarda, M.,
#' Blagotić, A., Shangguan, W., Wright, M. N., Geng, X., Bauer-Marschallinger, B.,
#' Guevara, M.A. Vargas, R., MacMillan, R. A., Batjes, N. H., Leenaars, J. G. B.,
#' Ribeiro, E., Wheeler, I., Mantel, S., & Kempen, B. (2017). SoilGrids250m:
#' Global gridded soil information based on machine learning. *PLoS one*, 12(2),
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
#'   \item{sand}{Sand fraction of soil by weight}
#'   \item{clay}{Clay fraction of soil by weight}
#'   \item{Dryprop}{Proportion of days that all layers within the MCS are dry when
#'    soil temperature at 50 cm >5 degrees Celsius}
#'   \item{CwetWinter}{Number of consecutive days with all MCS layers wet during
#'   the winter}
#'   \item{CdrySummer}{Number of consecutive days with all MCS layers dry during
#'   the summer}
#'   \item{Cwet8}{Number of consecutive days with any layer wet when soil
#'   temperature at 50 cm depth is >8 degrees Celsius}
#'   \item{Dryall}{Number of days with all MCS layers dry}
#'   \item{Dryany}{Number of days when any soil layer in the MCS is dry}
#' }
#' @source The cases in this data frame represent sites that were simulated for
#' Bradford et al. (2019). The simulations were conducted using mean current and
#' future climate conditions at 41,477 dryland sites across western North America
#' defined by a 10-km grid. The data for Wyoming represent 2126 of those cells.
#' For each 10-km cell, Bradford et al. (2019) simulated five different soil types:
#' site-specific soils derived from STATSGO (Miller and White 1998) and four fixed
#' soil types (see 'setsoiltypes' data for details), for a total of 10630 sites
#' with unique location x soil combinations in this dataset. The six output results
#' are from simulations under current climate conditions (1980-2010).
#'
#' Of the five soil types simulated for each 10 km cell, only the site-specific
#' soils varied with depth (the sand and clay content of the set soil types were
#' constant throughout the soil profile). We calculated depth-weighted sand and
#' clay content for each site-specific soil for the simulated sites.
#'
#' We calculated 30-year (1981 to 2010) normal monthly temperature and
#' precipitation data at 30-arcsecond resolution from DayMet (Thornton et al.,
#' 2018) using Google Earth Engine (Gorelick et al., 2017). We then calculated 19
#' bioclimatic variables using the \code{\link[dismo]{biovars}} function in R package dismo
#' (Hijmans, 2017). To get bioclim variables for the 10 km simulated sites, we
#' aggregated the 30-arcsecond bioclim variables to ~10 km and reprojected them
#' to Albers Equal Area.
#'
#' **References**\cr
#' Bradford, J. B., Schlaepfer, D. R., Lauenroth, W. K., Palmquist, K. A.,
#' Chambers, J. C., Maestas, J. D., & Campbell, S. B. (2019). Climate-driven
#' shifts in soil temperature and moisture regimes suggest opportunities to
#' enhance assessments of dryland resilience and resistance. *Front. Ecol. Evol.*
#' 7:358. https://doi.org/10.3389/fevo.2019.00358
#'
#' Gorelick, N., Hancher, M., Dixon, M.,  Ilyushchenko, S., Thau, D., & Moore, R.
#' (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone.
#' *Remote Sens. Environ.* 202, pp.18-27. https://doi.org/10.1016/j.rse.2017.06.031
#'
#' Hijmans, R. J., Phillips, S., Leathwick, J., & Elith, J. (2017). dismo: Species
#' Distribution Modeling. R package version 1.1-4. https://CRAN.R-project.org/package=dismo
#'
#' Miller, D. A., & White, R. A. (1998). A conterminous United States multilayer
#' soil characteristics dataset for regional climate and hydrology modeling.
#' *Earth Interact.* 2, pp.1–26.
#' https://doi.org/10.1175/1087-3562(1998)002<0001:ACUSMS>2.3.CO;2
#'
#' Schlaepfer, D. R., Lauenroth, W. K., & Bradford, J. B. (2012). Ecohydrological
#' niche of sagebrush ecosystems. *Ecohydrology* 5: 453-466.
#' https://doi.org/10.1002/eco.238
#'
#' Schlaepfer, D. R., & Andrews, C. A. (2018). rSFSW2: Simulation Framework for
#' SOILWAT2. R package version 3.0.0.
#'
#' Schlaepfer, D. R., & Murphy, R. (2018). rSOILWAT2: An Ecohydrological
#' Ecosystem-Scale Water Balance Simulation Model. R package version 2.3.2.
#'
#' Soil Survey Staff. (2014). Keys to Soil Taxonomy, 12th Edn. Washington, DC:
#' USDA-Natural Resources Conservation Service.
#'
#' Thornton, P. E., Thornton, M. M., Mayer, B. W., Wei, Y., Devarakonda, R., Vose,
#' R. S., & Cook, R. B. (2018). Daymet: Daily Surface Weather Data on a 1-km Grid
#' for North America, Version 3. ORNL DAAC, Oak Ridge, Tennessee, USA.
#' https://doi.org/10.3334/ORNLDAAC/1328
"subsetcells"





#' Four set soil types used for SOILWAT2 simulations
#'
#' A data frame of the percentage of sand and clay by weight of four set soil
#' types used for SOILWAT2 (Schlaepfer et al., 2012; Schlaepfer & Andrews, 2018;
#' Schlaepfer & Murphy, 2018) simulations in Bradford et al. (2019). Rownames
#' provide the treatment number for each soil type (these numbers are appended
#' as a decimal onto the unique site numbers 'site_id' in the 'subsetcells'
#' dataset to create the 'site_ids').
#'
#' @format A data frame with 4 rows and 2 variables:
#'
#' \describe{
#'   \item{sand}{Percentage sand content of soil by weight}
#'   \item{clay}{Percentage clay contet of soil by weight}
#' }
#'
#' @source Soil texture types were determined through an analysis of the soil
#' texture classes for big sagebrush ecosystems across the western US as described
#' in Appendix S3 of Palmquist et al. (2021).
#'
#' **References**\cr
#' Bradford, J. B., Schlaepfer, D. R., Lauenroth, W. K., Palmquist, K. A.,
#' Chambers, J. C., Maestas, J. D., & Campbell, S. B. (2019). Climate-driven
#' shifts in soil temperature and moisture regimes suggest opportunities to
#' enhance assessments of dryland resilience and resistance. *Front. Ecol. Evol.*
#' 7:358. https://doi.org/10.3389/fevo.2019.00358
#'
#' Palmquist, K. A., Schlaepfer, D. R., Renne, R. R., Torbit, S. C., Doherty, K. E.,
#' Remington, T. E., Watson, G., Bradford, J. B. and Lauenroth, W. K., 2021.
#' Divergent climate change effects on widespread dryland plant communities driven
#' by climatic and ecohydrological gradients. *Glob. Change Biol.*
#' https://doi.org/10.1111/gcb.15776
#'
#' Schlaepfer, D. R., Lauenroth, W. K., & Bradford, J. B. (2012). Ecohydrological
#' niche of sagebrush ecosystems. *Ecohydrology* 5: 453-466.
#' https://doi.org/10.1002/eco.238
#'
#' Schlaepfer, D. R., & Andrews, C. A. (2018). rSFSW2: Simulation Framework for
#' SOILWAT2. R package version 3.0.0.
#'
#' Schlaepfer, D. R., & Murphy, R. (2018). rSOILWAT2: An Ecohydrological
#' Ecosystem-Scale Water Balance Simulation Model. R package version 2.3.2.
"setsoiltypes"
