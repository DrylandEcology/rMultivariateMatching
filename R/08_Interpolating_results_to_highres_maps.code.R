##########################################
#
# Purpose: Interpolate simulation results to create high-resolution maps
#
# Inputs: "SubsetCells_bioclim+soils+results_WY_interpolation.csv"
#         "TargetCells_Bioclim+soils+matches_WY_interpolation.csv"
#         "drylands_wyoming.tif"
#
# Outputs: "Interpolated_output_vars_maps.png"
#
# Rachel Renne
# April 4, 2021
#
##########################################

################
# Step 1: Set up folders

# Folder where input and output data are stored
datafolder <- getwd()

# Folder to store figures
figurefolder <- getwd()

#########################
# Step 2: Read in data and interpolate results

# Read in bioclim and soils data for subset cells
subset <- read.csv(paste0(datafolder,"/SubsetCells_bioclim+soils+results_WY_interpolation.csv"), row.names = 1)
row.names(subset) <- paste0("00", row.names(subset))
# Just pull out results:
subset1x <- subset[,28:33]

# Read in data with matches
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

# Create dataframe of interpolated variables
resultsx <- data.frame(cbind(bioclim1[,1:2], subset1x[bioclim1$matches,]))

# now just include thoses with matching quality <= 1.5
resultsx1 <- resultsx[qual <= 1.5,]

# Read in Wyoming drylands raster
wydry <- raster(paste0(datafolder, "/drylands_wyoming.tif"))

# get Wyoming shapefile:
us <- getData("GADM", country = "USA", level = 1)
wy <- us[us$NAME_1 == "Wyoming",]
rm(us)

# Create dataframe of breaks for displaying differences in maps:
basebks <- list(0,0,0,c(0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1),
                c(0,15,25,35,45,60,75,90,105,122),
                c(23,30,35,40,45,60,75,90,105,122),
                c(0,15,25,35,45,70,95,120,145,164),
                c(33,80,120,160,200,240,280,320,360,366),
                c(46,60,70,80,90,135,180,225,270,315,366))


# Set up variable names:
vars <- c("","","","A) Dryprop", "B) CwetWinter","C) CdrySummer", "D) Cwet8", "E) Dryall", "F) Dryany")

# Set up colors
# For cool/wet variables
cols <- c(rev(c("#fecc5c","#fd8d3c","#f03b20","#bd0026")),
          c("#bdd7e7","#6baed6","#3182bd","#08519c"))

# cool/wet with one extra 0n cool/wet side
colsx <- c(rev(c("#fecc5c","#fd8d3c","#f03b20","#bd0026")),
           c("#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))

# for cdrysummer
colsxy <- c(rev(c("#bdd7e7","#6baed6","#3182bd","#08519c")),
            c("#ffeda0","#fed976","#feb24c","#bd0026","#800026"))

# For hot/dry variables?
cols1 <- rev(c(rev(c("#fecc5c","#fd8d3c","#f03b20","#bd0026")),
               c("#bdd7e7","#6baed6","#3182bd","#08519c")))

# For Dry all
dall <- c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","#000000")

# dryany
dany <- rev(c(rev(c("#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026")),
              c("#bdd7e7","#6baed6","#3182bd","#08519c")))

# Set up cols:
colx <- list(0,0,0,cols1, colsx, colsxy, colsx, dall, dany)

par(mfrow = c(2,3))
for (i in 4:9){
  plot(1:length(colx[[i]])~1, col = colx[[i]], pch = 16, cex = 2)
}

# Set up png file:
png(file = paste0(figurefolder, "/Interpolated_output_vars_maps.png"), width = 18, height = 12, units = 'in', res = 600)
par(mfrow = c(2,3), mar = c(1,1,2,1), mgp = c(3,0.7,0))
# Rasterize differences & plot:
for (i in 4:9){
  print(paste0("Now working on ", vars[i],"."))
  # Create spatial points dataframe from differences:
  ptsx <- SpatialPointsDataFrame(resultsx1[,1:2], data = data.frame(resultsx1[,i-1]), proj4string = crs(wydry))
  # Rasterize
  r <- rasterize(ptsx, wydry, field = names(ptsx), fun = mean, na.rm = T)
  
  # Calculate breaks
  bks <- basebks[[i]]
  truebks <- basebks[[i]]
  truecols <- colx[[i]]
  
  # Set up axis label locations
  axisat <- -111
  for (ii in 2:(length(bks))){
    axisat <- c(axisat, (axisat[ii-1]+(6.9/length(truecols))))
  }
  
  # create figures
  image(r, col = truecols, breaks = bks, bty = "n", xaxt = "n",yaxt="n",ylim = c(40.325,45.01),
        main = "",
        xlab = "", ylab = "")
  plot(wy, add = T, lwd = 2)
  for (xx in 1:length(truecols)){
    polygon(x = c(axisat[xx],axisat[xx],axisat[xx+1],axisat[xx+1]), y = c(40.75,40.95,40.95,40.75),border = truecols[xx], col = truecols[xx])
  }
  polygon(x = c(-111, -111, -104.1 ,-104.1), y = c(40.75,40.95,40.95,40.75), lwd = 1.5)
  mtext(vars[i], side = 3, line = 0, cex = 2, adj = 0)
  par(tcl = -0.3)
  axis(side = 1, line = -3.5, at = seq(-111,-104.1, length.out = length(bks)),cex.axis = 2,
       labels = bks)
  
}
#dev.off()
