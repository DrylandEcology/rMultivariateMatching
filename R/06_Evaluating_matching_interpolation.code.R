##########################################
#
# Purpose: Evaluate matching for interpolation case study
#
# Summary: Three methods to evaluate matching for site selection method
#           1) Matching quality variable: weighted Euclidean distance between Subset and Target cells
#           2) Standard deviation of differences between Subset and Target cells for a set of relevant variables
#           3) Geographic distance between Target and matched Subset cells and average geographic distance
#              between Subset cells matched to a Target cells and the eight adjacent matched Subset cells
#
# Inputs: "drylands_wyoming.tif"
#         "TargetCells_Bioclim+soils+matches_WY_interpolation.csv"
#
# Outputs: "Mapped Coverage_interpolation.png"
#          "Mapped Coverage_soils_interpolation.png"
#          "Variable Coverage_interpolation.png"
#          "SD_Diffs_allcells+matched_interpolation.png"
#          "Distances_subset+target_cells_interpolation.png"
#          "Distances_adjacent_subset_cells_interpolation.png"
#
# Rachel Renne
# March 27, 2021
#
##########################################

library(raster)
library(distances)

#####################################################################################
# Step 1: Calculate matching quality (weighted Euclidean distance between target and matched subset cells)

# Folder where input and output data are stored
datafolder <- getwd()

# Folder to store figures
figurefolder <- getwd()

# Read in bioclim and soils with matches data
bioclim <- readRDS(paste0(datafolder,"/TargetCells_Bioclim+soils+matches_WY_interpolation.csv"))
# Note: cellnumbers in this dataframe correspond to the unique cellnumbers associated with each
# target cell. These serve as unique identifiers throughout the code and can
# be extracted from the raster datasets that define the target cells using the 
# extract function in the raster package with "cellnumbers" = TRUE

# limit to just the 6 variables of interest (i.e., 6 bioclim matching variables):
bioclim1a <- data.frame(bioclim[,c(1,2,3,6,11,14,17,20,22:24)])

# transform those variables by dividing by matching criterion to make 
# Euclidean distance of "1" for each variable equal to the maximum acceptable 
# difference between target and subset cells

# Matching criteria in parentheses
# "bioclim_01" = mean annual temperature (0.7)
# "bioclim_04" = temperature seasonality (42)
# "bioclim_09" = mean temperature of the driest quarter (3.3)
# "bioclim_12" = mean annual precipitation (66)
# "bioclim_15" = precipitation seasonality (5.4)
# "bioclim_18" = precipitation of the warmest quarter (18.4)
# "clay" = depth-weighted percentage clay (0.053)
# "sand" = depth-weighted percentage sand (0.038)

# Standardize variables of interest according to percentage of range
bioclim1 <- cbind(bioclim1a[,1:2],bioclim1a$bioclim_01/0.7, 
                  bioclim1a$bioclim_04/42, 
                  bioclim1a$bioclim_09/3.3,
                  bioclim1a$bioclim_12/66,
                  bioclim1a$bioclim_15/5.4, 
                  bioclim1a$bioclim_18/18.4,
                  bioclim1a$sand/0.053,
                  bioclim1a$clay/0.038) 
# fix column names
colnames(bioclim1) <- c(colnames(bioclim1a)[1:2],"bioclim_01","bioclim_04","bioclim_09", "bioclim_12","bioclim_15","bioclim_18","sand","clay")

# Add matches onto bioclim
bioclim1$matches <- paste0("00", bioclim1a$matches)

# Read in Subset cells with full climate x soil data
subset <- read.csv(paste0(datafolder,"/SubsetCells_bioclim+soils+results_WY_interpolation.csv"), row.names = 1)
row.names(subset) <- paste0("00", row.names(subset))

# create df of variables of interest for subset cells
subset1a <- subset[,c(2,3,4,7,10,15,18,21,24,26,27)]

# Standardize variables of interest according to percentage of range
subset1 <- cbind(subset1a[,1:2],subset1a$bioclim_01/0.7, 
                  subset1a$bioclim_04/42, 
                  subset1a$bioclim_09/3.3,
                  subset1a$bioclim_12/66,
                  subset1a$bioclim_15/5.4, 
                  subset1a$bioclim_18/18.4,
                  subset1a$sand/0.053,
                  subset1a$clay/0.038) 
# fix column names
colnames(subset1) <- c(colnames(subset1a)[1:2],"bioclim_01","bioclim_04","bioclim_09", "bioclim_12","bioclim_15","bioclim_18","sand","clay")

# Use matches column in bioclim1 to calculate weighted Euclidean distance between target and matched subset cells
d1 <- (subset1[bioclim1$matches,3]-bioclim1[,3])^2
for (cv in 2:(ncol(bioclim1)-3)){
  d1<- cbind(d1,(subset1[bioclim1$matches,2+cv]-bioclim1[,2+cv])^2)
}
sum6 <- apply(d1[,1:6],1, sum)

# Final weighted Euclidean distance between each Target cell and its matched Subset cell (i.e. matching quality variable):  
qual <- data.frame(distance = sqrt(sum6))

# Create soil matching quality variables
soilqual <- data.frame(distance = sqrt(apply(d1[,7:8],1,sum)))

############################
# Map matching quality

# Read in Wyoming drylands raster
wydry <- raster(paste0(datafolder,"/drylands_wyoming.tif"))

# Create spatial points dataframe
ptsx <- SpatialPointsDataFrame(bioclim1[,1:2], data = qual, proj4string = crs(wydry))

# Create raster of qual (matching quality) using wydry as a template
r <- rasterize(ptsx, wydry, field = qual$distance, fun = mean)

# get Wyoming shapefile:
us <- getData("GADM", country = "USA", level = 1)
wy <- us[us$NAME_1 == "Wyoming",]
rm(us)

# Designate colors, breaks, and plotting settings
cols = rev(c("#d7191c","#fdae61","#abd9e9","#2c7bb6"))
bks = c(0,0.5,1,1.5,13)
par(mar = c(1,1,2,1), tcl = -0.5, mgp = c(3,1,0))

#png(file = paste0(figurefolder, "/Mapped Coverage_interpolation.png"), width = 6, height = 6, units = 'in', res = 300)
image(r, col = cols, breaks = bks, bty = "n", xaxt = "n",yaxt="n",ylim = c(40.325,45.01),
      main = "",
      xlab = "", ylab = "")
plot(wy, add = T, lwd = 2)
points(subset1a$Y_WGS84 ~ subset1a$X_WGS84, pch = 16, cex = 0.2, col = "black")
mtext("Matching quality", 3, cex = 1.6)
legend("bottom" , legend = c("0 to 0.5","0.5 to 1","1 to 1.5","1.5 to 12.4"), fill = cols,
       bty = "n", cex = 1, ncol = 4, x.intersp = 0.6)
#dev.off()

############
# Map soil matching quality variable:
# Create spatial points dataframe
ptsx <- SpatialPointsDataFrame(bioclim1[,1:2], data = soilqual, proj4string = crs(wydry))

# Create raster of qual (matching quality) using wydry as a template
r <- rasterize(ptsx, wydry, field = soilqual$distance, fun = mean)

# Set up breaks
bks = c(0,0.5,1,1.5,13)

#png(file = paste0(figurefolder, "/Mapped Coverage_soils_interpolation.png"), width = 6, height = 6, units = 'in', res = 300)
image(r, col = cols, breaks = bks, bty = "n", xaxt = "n",yaxt="n",ylim = c(40.325,45.01),
      main = "",
      xlab = "", ylab = "")
points(subset1a$Y_WGS84 ~ subset1a$X_WGS84, pch = 16, cex = 0.2, col = "black")
plot(wy, add = T, lwd = 2)
mtext("Matching quality", 3, cex = 1.6)
legend("bottom" , legend = c("0 to 0.5","0.5 to 1","1 to 1.5","1.5 to 4.9"), fill = cols,
       bty = "n", cex = 1, ncol = 4, x.intersp = 0.6)
#dev.off()

########################
# Look at coverage of matching variables (scatterplots)

# Scatterplot of variables:
# bioclim_01 = Mean Annual Temperature
# bioclim_04 = Temperature Seasonality
# bioclim_09 = Mean Temp of Driest Quarter
# bioclim_12 = Mean Annual Precipitation
# bioclim_15 = Precipitation seasonality
# bioclim_18 = Precip of warmest quarter

#png(file = paste0(figurefolder,"/Variable Coverage_interpolation.png"),width = 9, height = 3, units = 'in', res = 300)
par(mar = c(3,3,1,1), mgp = c(1.5,0.3,0), tcl = 0.5, mfrow = c(1,3))
# MAP~MAT
smoothScatter(bioclim1a$bioclim_12~bioclim1a$bioclim_01, xlab = "Mean annual temperature (C)",
              ylab = "Mean annual precipitation (mm)", 
              colramp = colorRampPalette(c("white","#a1dab4","#41b6c4","#2c7fb8","#253494")),
              nrpoints = 0, nbin = 300, cex.lab = 1.45, cex.axis = 1.5)
points(subset1a$bioclim_12 ~ subset1a$bioclim_01, pch = 16, cex = 0.5)

# Pseas~Tseas
smoothScatter(bioclim1a$bioclim_15~bioclim1a$bioclim_04, xlab = "Temperature seasonality",
              ylab = "Precipitation seasonality", 
              colramp = colorRampPalette(c("white","#a1dab4","#41b6c4","#2c7fb8","#253494")),
              nrpoints = 0, nbin = 300, cex.lab = 1.45, cex.axis = 1.5)
points(subset1a$bioclim_15~subset1a$bioclim_04, pch = 16, cex = 0.5)

# Pwarmq~TdryQ
smoothScatter(bioclim1a$bioclim_09~bioclim1a$bioclim_18, xlab = "Precipitation of warmest quarter (mm)",
              ylab = "Temperature of driest quarter (C)", 
              colramp = colorRampPalette(c("white","#a1dab4","#41b6c4","#2c7fb8","#253494")),
              nrpoints = 0, nbin = 300, cex.lab = 1.45, cex.axis = 1.5)
points(subset1a$bioclim_09~subset1a$bioclim_18, pch = 16, cex = 0.5)
#dev.off()


#####################################################################################
# Step 2: Calculate standard deviation of differences between subset and target cells
# for a set of cells relevant to the project

# Make bioclim dataset that excludes sites that are have matching quality variable >1.5 
#(weighted Euclidean distance):
marea <- bioclim[qual <= 1.5, ]

# Create results dataframe to store standard deviations:
results <- data.frame(area = c('allcells','matchedcells'),
                      a = NA, b = NA, c = NA, d = NA, e = NA, f = NA, g = NA, h = NA, i = NA,
                      j = NA, k = NA, l = NA, m = NA, n = NA, o = NA, p = NA, q = NA, r = NA, s = NA,
                      t = NA, u = NA)
colnames(results)[2:22] <- colnames(bioclim)[3:23]

# Calculate SD of diffs for all and for matched cells:
for (i in 3:23){
  results[1,i-1] <- sd(subset[bioclim1$matches,i+4]-bioclim[,i])
  results[2,i-1] <- sd(subset[bioclim1[qual$distance < 1.5,]$matches,i+4]-marea[,i])
}

# Convert to percentage sand & clay
results$sand <- results$sand*100
results$clay <- results$clay*100

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

# Create barplot showing standard deviation of differences between Target and Subset cells
# For all cells and for only those cells with matching quality <= 1.5 ("matched")
#png(paste0(figurefolder,"/SD_Diffs_allcells+matched_interpolation.png"), width = 5, height = 7, units = "in", res = 300)
par(mfrow = c(1,1))
par(mar = c(3,10,3,1), mgp = c(3,0.3,0), tcl = 0.3, lwd =1)
barplot(as.matrix(results[,22:2]), beside = T, horiz = T, col = rep(c("grey",0),19),
        names.arg = rev(bioclim_vars), las = 1, cex.names = 0.7,
        main = "Standard deviation of differences", xlim = c(0,62))
legend(x = 15, y = 26, legend = c("Matched only","All cells"), fill= c(0,"grey"), bty = "n")
box()
#dev.off()

#####################################################################################
# Step 3: Calculate distance between target and matched subset cells 
# and distance between subset cells matched to each target cells and the subset cells matched to the
# eight adjacent neighbors

####################
# Distance between Target and matched Subset cells

# Create matrix of lat/long for cells with matching quality <= 1.5
# First two columns are coordinates of matched Subset cells
# Second two columns are coordinates of Target cells
pts1 <- cbind(subset[paste0("00",marea[,24]),3:4], marea[,1:2])

# Calculate distances between pairs of points
alldistsqual <- pointDistance(p1 = pts1[,1:2], p2 = pts1[,3:4],lonlat = T)

# Convert to km
distkm <- alldistsqual/1000
summary(distkm)
hist(distkm, breaks = 100)

###################
# Create map of distances:
distmatrix <- cbind(pts1[,3:4], distkm)
rownames(distmatrix) <- c(1:nrow(distmatrix))
colnames(distmatrix)[3] <- "distance"

# Create spatial points dataframe from differences:
ptsx <- SpatialPointsDataFrame(distmatrix[,1:2], data = data.frame(dist = distkm), proj4string = crs(wydry))

# Rasterize using wydry as a template
r <- rasterize(ptsx, wydry, field = ptsx$dist, fun = mean)

# Designate colors
truecols <- c("#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#993404","#662506")

# set breaks
bks <- c(0, 6, 12, 18, 36, 72, 144, 288, round(max(distkm)))

#png(file = paste0(figurefolder,"/Distances_subset+target_cells_interpolation.png"), width = 6, height = 6, units = 'in', res = 300)
par(mfrow = c(1,1), mar = c(1,1,1,1), mgp = c(3,0.25,0))

# Set up axis label locations
axisat <- -111
for (ii in 2:(length(bks))){
  axisat <- c(axisat, (axisat[ii-1]+(6.9/length(truecols))))
}

# create figures
image(r, 
      col = truecols, 
      breaks = bks, 
      ylim = c(40.3, 45.01),
      xlim = c(-111.06, -104.04), useRaster = T,
      xlab = "", ylab ="",
      bty = "n", xaxt = "n",yaxt="n")
plot(wy, add = T, lwd = 2)
for (xx in 1:length(truecols)){
  polygon(x = c(axisat[xx],axisat[xx],axisat[xx+1],axisat[xx+1]), y = c(40.75,40.95,40.95,40.75),border = truecols[xx], col = truecols[xx])
}
polygon(x = c(-111, -111, -104.1 ,-104.1), y = c(40.75,40.95,40.95,40.75), lwd = 1.5)
mtext("Distance to Subset cell (km)", side = 1, line = -0.6, cex = 1)
par(tcl = -0.3)
axis(side = 1, line = -2.5, at = seq(-111,-104.1, length.out = length(bks)),cex.axis = 0.9,
     labels = bks)
#dev.off()  

##################################################
# Calculate distance between subset cell matched to each target cells 
# and the subset cells matched to its 8 adjacent neighbors

# Make points matrix where cols 3&4 are coordinates of subset cell and 1&2 of target cell
pts3 <- cbind(marea[,1:2], subset[paste0("00",marea[,24]),3:4])

# Create points dataframe
ptsx <- SpatialPointsDataFrame(pts3[,1:2], data = data.frame(thispt = rep(1,nrow(pts3))), proj4string = crs(wydry))

# Rasterize using wydry as template
r <- rasterize(ptsx, wydry, field = 1 , fun = 'first')

# Use adjacent function to get up to 8 adjacent cells
# First column gives the focal cell, second column gives any non-NA cell it is adjacent to.
a <- adjacent(r, cells=as.numeric(row.names(ptsx)), directions=8, pairs=TRUE)

# Remove nonexistant neighbors (i.e., neighbors not included in marea):
a1 <- a[which(as.character(a[,2]) %in% row.names(ptsx)),]

# Calculate distances between subset cells matched to each set of neighbors
dista <- pointDistance(p1 = pts3[as.character(a1[,1]),3:4], p2 = pts3[as.character(a1[,2]),3:4],lonlat = T)

# Calculate average distance between subset cell matched to target cell and subset cells matched to its 8 neighbors:
avgdist <- tapply(dista, as.factor(a1[,1]), FUN = mean)

# Look for cells without neighbors (just FYI):
nrow(pts3) - length(avgdist)

# Create distances dataframe:
adjdist <- data.frame(cellnumber = names(avgdist), avgdist_km = avgdist/1000)
# Add in coordinates
adjdist$x <- pts3[adjdist$cellnumber, 1]
adjdist$y <- pts3[adjdist$cellnumber, 2]

################################
# Map these distances:

# Create spatial points dataframe from differences:
ptsx <- SpatialPointsDataFrame(adjdist[,3:4], data = data.frame(dist = adjdist[,2]), proj4string = crs(wydry))

# Rasterize using wydry as template
r <- rasterize(ptsx, wydry, field = ptsx$dist, fun = mean)

# Designate colors
truecols <- c("#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#993404","#662506")

# Designate breaks
bks <- c(0, 6, 12, 18, 36, 72, 144, 288, round(max(adjdist$avgdist_km)))

# Create file
#png(file = paste0(figurefolder,"/Distances_adjacent_subset_cells_interpolation.png"), width = 6, height = 6, units = 'in', res = 300)
par(mfrow = c(1,1), mar = c(1,1,1,1), mgp = c(3,0.25,0))

# Set up axis label locations
axisat <- -111
for (ii in 2:(length(bks))){
  axisat <- c(axisat, (axisat[ii-1]+(6.9/length(truecols))))
}

# create figures
image(r, 
      col = truecols, 
      breaks = bks, 
      ylim = c(40.3, 45.01),
      xlim = c(-111.06, -104.04), useRaster = T,
      xlab = "", ylab ="",
      bty = "n", xaxt = "n",yaxt="n")
plot(wy, add = T, lwd = 2)
for (xx in 1:length(truecols)){
  polygon(x = c(axisat[xx],axisat[xx],axisat[xx+1],axisat[xx+1]), y = c(40.75,40.95,40.95,40.75),border = truecols[xx], col = truecols[xx])
}
polygon(x = c(-111, -111, -104.1 ,-104.1), y = c(40.75,40.95,40.95,40.75), lwd = 1.5)
mtext("Distance (km)", side = 1, line = -0.6, cex = 1)
par(tcl = -0.3)
axis(side = 1, line = -2.5, at = seq(-111,-104.1, length.out = length(bks)),cex.axis = 0.9,
     labels = bks)
#dev.off()
