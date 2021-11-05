#' Plot a raster with a custom legend below
#'
#' Plots a raster with binned colors and adds a legend below the image.
#'
#'
#' @param x raster. A raster to be plotted.
#'
#' @param thisVariable character. The name of the variable in the raster. Defaults
#' to the name of the raster.
#'
#' @param round_dec numeric. The number of decimal places to round the labels
#' of the legend. Defaults to 0.
#'
#' @param cols vector of colors to use in the legend. Default colors range from
#' yellow to dark brown through 8 distinct colors.
#'
#' @param bks vector of breaks to use in designating which values get assigned to
#' each color in `cols`. Must have one more element than `cols`. Unless `bks` are
#' designated, this vector will be calculated internally from `cols`.
#'
#' @return a plot of the raster with a legend.
#'
#'
#' @examples
#' # Load bioclim data for Target Cells (from rMultivariateMatchingAlgorithms package)
#' data(bioclim)
#' legendPlot(bioclim[[1]])
#'

loocv <- function(matches){

#####################################################################################
# Step 1: Calculate matching quality (weighted Euclidean distance between target and matched subset cells)

# Folder where input and output data are stored
datafolder <- getwd()

# Folder to store figures
figurefolder <- getwd()

# Read in bioclim and soils data for subset cells
subset <- read.csv(paste0(datafolder,"/SubsetCells_bioclim+soils+results_WY_interpolation.csv"), row.names = 1)
row.names(subset) <- paste0("00", row.names(subset))

#pull out bioclim+soils data:
bioclimx <- subset[,c(3,4,7:27)]

# Pull out just sites with site-specific soils (every 5th site)
bioclim1 <- bioclimx[seq(1,nrow(bioclimx), by = 5),]

# Transform
bioclim2 <- cbind(bioclim1[,1:2],bioclim1$bioclim_01/0.7,
                  bioclim1$bioclim_04/42,
                  bioclim1$bioclim_09/3.3,
                  bioclim1$bioclim_12/66,
                  bioclim1$bioclim_15/5.4,
                  bioclim1$bioclim_18/18.4,
                  bioclim1$sand/0.053,
                  bioclim1$clay/0.038)
colnames(bioclim2) <- c("x","y","MAT","MAP","ppt.seas","temp.seas","tdryq","pwarmq","sand","clay")

# Calculate distance matrix for all cells:
xdist <- distances(bioclim2[,3:8], id_variable = row.names(bioclim2))
cell_numbers <- rownames(bioclim2)

# Find 2 nearest neighbors (first will be self) using Euclidean distance of weighted climate variables
neighbors <- nearest_neighbor_search(xdist, k = 2, search_indices = 1:2126,
                                     query_indices = 1:2126)
neighbors1 <- cell_numbers[t(neighbors[2,])]

####################
# Step 3: Find best match among available soils for each climatic match

# Calculate the distance between site-specific soils and each matched cell:
# Get matched values:
site_vars1a <- bioclim1[t(neighbors[2,]),]

# Calculate distance between site specific soils and matched cells
site_dist <- sqrt((site_vars1a[,22]/0.053-bioclim1[,22]/0.053)^2+(site_vars1a[,23]/0.038-bioclim1[,23]/0.038)^2)
# Set names to the Subset cell name that corresponds to site-specific soils
names(site_dist) <- neighbors1


# Create dataframe of set soil type
sand <- c(0.3, 0.27, 0.66, 0.16)
clay <- c(0.18, 0.35, 0.09, 0.09)
soils <- cbind(sand, clay)
rownames(soils) <- c("2","3","4","5")

# Combine matched cells soils and allcells soils:
allsoils <- rbind(soils,bioclim1[,22:23])

# Transform soils:
allsoils1 <- cbind(allsoils[,1]/0.053, allsoils[,2]/0.038)
rownames(allsoils1) <- c("2","3","4","5",row.names(bioclim1))

# Calculate distances between set soil types and matched cells:
xdist1 <- distances(allsoils1, id_variable = row.names(allsoils1))
cell_numbers <- rownames(allsoils1)
neighbors_soils <- nearest_neighbor_search(xdist1, k = 1, search_indices = c(1:4),
                                           query_indices = c(5:nrow(allsoils1)))

# Create dataframe with matched simiulated site variables corresponding to each cell
site_vars2a <- allsoils1[neighbors_soils,]

# Calculate distance between site specific soils and matched cells
setsoil_dist <- sqrt((site_vars2a[,1]-bioclim1[,22]/0.053)^2+(site_vars2a[,2]-bioclim1[,23]/0.038)^2)

# Determine min.dist between site_specific distance and min dist to set soil types:
min_dist <- apply(cbind(site_dist, setsoil_dist), 1, FUN = which.min)

# Create vector of rownames from setsoil_dist that reflects unique ids:
setsoil_names <- paste0("00",as.integer(as.numeric(names(site_dist))), ".",names(setsoil_dist))
names_matrix <- cbind(names(site_dist), setsoil_names)
names_selection_matrix <- cbind(1:length(min_dist), min_dist)
matches <- names_matrix[names_selection_matrix]

# Create a matches column for target cells
bioclimx1 <- data.frame(cbind(bioclim1,matches))

# Add matches on for set soils:
clim_match <- as.integer(as.numeric(matches))

# Pull out each soil type:
bioclimx2 <- bioclimx[seq(2,nrow(bioclimx), by = 5),]
bioclimx2$matches <- paste0("00", clim_match,".2")

bioclimx3 <- bioclimx[seq(3,nrow(bioclimx), by = 5),]
bioclimx3$matches <- paste0("00", clim_match,".3")

bioclimx4 <- bioclimx[seq(4,nrow(bioclimx), by = 5),]
bioclimx4$matches <- paste0("00", clim_match,".4")

bioclimx5 <- bioclimx[seq(5,nrow(bioclimx), by = 5),]
bioclimx5$matches <- paste0("00", clim_match,".5")

# Put all bioclim+soils data together with matches:
alldata <- rbind(bioclimx1, bioclimx2, bioclimx3, bioclimx4, bioclimx5)

##########################################
# Step 4: Calculate standard deviation of diffs for clim/soil variables

# Read in list of variable names
bioclim_vars <- c("Annual Mean Temperature",
                  "Mean Diurnal Range",
                  "Isothermality (BIO2/BIO7)",
                  "Temperature Seasonality (SD*100)",
                  "Max Temp of Warmest Month",
                  "Min Temp of Coldest Month",
                  "Temp Annual Range",
                  "Mean Temp of Wettest Quarter",
                  "Mean Temp of Driest Quarter",
                  "Mean Temp of Warmest Quarter",
                  "Mean Temp of Coldest Quarter",
                  "Annual Precipitation",
                  "Precip of Wettest Month",
                  "Precip of Driest Month",
                  "Precip Seasonality (Coef of Var)",
                  "Precip of Wettest Quarter",
                  "Precip of Driest Quarter",
                  "Precip of Warmest Quarter",
                  "Precip of Coldest Quarter",
                  "Sand fraction","Clay fraction")

# Results df for sd of differences:
results_sd <- data.frame(variable = bioclim_vars, six_vars = NA)

# Record sd and max diffs in results
for (i in 3:23){
  results_sd[i-2,2] <- sd(alldata[alldata$matches,i]-alldata[,i])
}

# Convert soils to percentage so it displays better among other variables
results_sd1 <- results_sd
results_sd1[20:21,2] <- c(results_sd[20:21,2]*100 )

# Visualize differences
# Visualize results:
#png(paste0(figurefolder,"/SD_Diffs_LOOCV_interpolation.png"), width = 5, height = 7, units = "in", res = 300)
par(mfrow = c(1,1))
par(mar = c(3,10,3,1), mgp = c(3,0.3,0), tcl = 0.3, lwd =1)
barplot(rev(c(t(results_sd1[,2]))), beside = T, horiz = T,
        names.arg = rev(c(bioclim_vars[1:19],"Sand (%)","Clay (%)")), las = 1, cex.names = 0.7,
        main = "Standard deviation of differences", xlim = c(0,13))
box()
#dev.off()

##########################################
# Step 5: Estimate matching errors

# Pull results in from subset dataframe & reorder to match alldata
results1 <- subset[,28:33]
results <- results1[row.names(alldata),]
# add on matches
results$matches <- alldata$matches


################################################################################

titles <- list(parse(text = paste0('"A)" ', ' ~ DRY[PROP]')), parse(text = paste0('"B)" ', ' ~ CWET[WINTER]')),
               parse(text = paste0('"C)" ', ' ~ CDRY[SUMMER]')), parse(text = paste0('"D)" ', ' ~ CWET[8]')),
               parse(text = paste0('"E)" ', ' ~ DRY[ALL]')), parse(text = paste0('"F)" ', ' ~ DRY[ANY]')))
rd <- c(2,0,0,0,0,0)
rd2 <- c(3,0,0,0,0,0)
rd3 <- c(3,1,1,1,1,1)

#png(paste0(figurefolder,"/Hist_Estimated_Errors_LOOCV_interpolation.png"), width = 9, height = 6, units = "in", res = 300)

par(mar = c(3,3,2,1), mfrow = c(2,3), mgp = c(1.7,0.3,0), tcl = 0.1)

for (i in 1:6){
  thisvar = round(sum((-results[results$matches,i]+results[,i])^2)/nrow(results),rd2[i])
  thismin = round(min(-results[results$matches,i]+results[,i]),rd[i])
  thismax = round(max(-results[results$matches,i]+results[,i]),rd[i])
  hist(-results[results$matches,i]+results[,i], breaks = 100,
       main = "",
       xlab = "Simulated - Matched value", cex.lab = 1.3, cex.axis = 1.2, ylim = c(0,1550))
  legend("topright", legend = c(paste0("min = ", thismin), paste0("max = ", thismax),
                                paste0("mean = ", round(mean(-results[results$matches,i]+results[,i],na.rm = T),rd3[i])),
                                parse(text = paste0('', ' ~ CV[error] == ', thisvar))), bty = "n", cex = 1.2)
  box()
  abline(h = 0)
  mtext(titles[[i]], side = 3, line = 0, adj = 0, cex = 1.2)
}
#dev.off()
}
