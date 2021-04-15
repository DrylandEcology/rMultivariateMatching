##########################################
#
# Purpose: Match Target cells to Subset cells simulated on a 10 km grid
#
# Summary: Each Subset cell was simulated with 5 different soils, and matching was best achieved
#          by first matching Target cells to Subset cells using climate, then choosing the best
#          match among the 5 simulated soil types
#
# Inputs: "TargetCells_Bioclim+soils_WY.csv"
#         "SubsetCells_bioclim+soils+results_WY_interpolation.csv"
#
# Outputs: "TargetCells_Bioclim+soils+matches_WY_interpolation.csv"
#
# Rachel Renne
# March 24, 2021
#
##########################################

library(distances)

################
# Step 1: Set up folders

# Folder where input and output data are stored
datafolder <- getwd()

##############################
# Step 2: Find best match using climate data

# Read in bioclim and soils data
bioclim <- readRDS(paste0(datafolder,"/TargetCells_Bioclim+soils_WY.csv"))
# Note: cellnumbers in this dataframe correspond to the unique cellnumbers associated with each
# target cell. These serve as unique identifiers throughout the code and can
# be extracted from the raster datasets that define the target cells using the 
# extract function in the raster package with "cellnumbers" = TRUE

#reverse soils data and remove cellnumbers column
bioclimx <- bioclim[,c(1:22,24,23)][,-1]

# Transform percentage sand and clay to fraction
bioclimx[,22] <- bioclimx[,22]/100
bioclimx[,23] <- bioclimx[,23]/100


# Calculate 10% of the range of sand and clay fraction for matching criteria
(max(bioclimx[,22])-min(bioclimx[,22]))*0.1
# 0.053
(max(bioclimx[,23])-min(bioclimx[,23]))*0.1
# 0.038

# Read in Subset cells
sites <- read.csv(paste0(datafolder,"/SubsetCells_bioclim+soils+results_WY_interpolation.csv"))

# Remove unneeded coordinates and other columns
sites1 <- sites[,c(4,5,8:28)]
colnames(sites1)[1:2] <- c("x","y")

# Create unique rownames for sites
row.names(sites1) <- paste0("00",sites$site_ids)

# Pull out just one set of sites to avoid climate duplicates
# Take every 5th site to have just those with site-specific soils
sites1a <- sites1[seq(1,nrow(sites1), by = 5),]

# Combine Subset and Target cells
bioclim1a <- rbind(sites1a, bioclimx)

bioclim2 <- cbind(bioclim1a[,1:2],bioclim1a$bioclim_01/0.7, 
                  bioclim1a$bioclim_04/42, 
                  bioclim1a$bioclim_09/3.3,
                  bioclim1a$bioclim_12/66,
                  bioclim1a$bioclim_15/5.4, 
                  bioclim1a$bioclim_18/18.4,
                  bioclim1a$sand/0.053,
                  bioclim1a$clay/0.038)
colnames(bioclim2) <- c("x","y","MAT","MAP","ppt.seas","temp.seas","tdryq","pwarmq","sand","clay")

# Calculate distance matrix for all cells:
xdist <- distances(bioclim2[,3:8], id_variable = row.names(bioclim2))
cell_numbers <- rownames(bioclim2)

# Find nearest neighbor using Euclidean distance of weighted climate variables
neighbors <- nearest_neighbor_search(xdist, k = 1, search_indices = 1:2126,
                                     query_indices = 2127:nrow(bioclim2))
neighbors1 <- cell_numbers[t(neighbors)]

####################
# Step 3: Find best match among available soils for each climatic match

# Calculate the distance between site-specific soils and each matched cell:
# Get matched values:
site_vars1a <- sites1a[t(neighbors),]

# Calculate distance between site specific soils and matched cells
site_dist <- sqrt((site_vars1a[,22]/0.053-bioclimx[,22]/0.053)^2+(site_vars1a[,23]/0.038-bioclimx[,23]/0.038)^2)
# Set names to the Subset cell name that corresponds to site-specific soils
names(site_dist) <- neighbors1


# Create dataframe of set soil type
sand <- c(0.3, 0.27, 0.66, 0.16)
clay <- c(0.18, 0.35, 0.09, 0.09)
soils <- cbind(sand, clay)
rownames(soils) <- c("2","3","4","5")

# Combine matched cells soils and allcells soils:
allsoils <- rbind(soils,bioclimx[,22:23])

# Transform soils:
allsoils1 <- cbind(allsoils[,1]/0.053, allsoils[,2]/0.038)

# Calculate distances between set soil types and matched cells:
xdist1 <- distances(allsoils1, id_variable = row.names(allsoils1))
cell_numbers <- rownames(allsoils1)
neighbors_soils <- nearest_neighbor_search(xdist1, k = 1, search_indices = c(1:4),
                                           query_indices = c(5:nrow(allsoils1)))

# Create dataframe with matched simiulated site variables corresponding to each cell
site_vars2a <- allsoils1[neighbors_soils,]

# Calculate distance between site specific soils and matched cells
setsoil_dist <- sqrt((site_vars2a[,1]-bioclimx[,22]/0.053)^2+(site_vars2a[,2]-bioclimx[,23]/0.038)^2)

# Determine min.dist between site_specific distance and min dist to set soil types:
min_dist <- apply(cbind(site_dist, setsoil_dist), 1, FUN = which.min)

# Create vector of rownames from setsoil_dist that reflects unique ids:
setsoil_names <- paste0("00",as.integer(as.numeric(names(site_dist))), ".",names(setsoil_dist))
names_matrix <- cbind(names(site_dist), setsoil_names)
names_selection_matrix <- cbind(1:length(min_dist), min_dist)
matches <- names_matrix[names_selection_matrix]

# Create a matches column for target cells
bioclimx1 <- cbind(bioclimx,as.numeric(matches))
colnames(bioclimx1)[24] <- "matches"

# Save file for next steps:
#saveRDS(bioclimx1, paste0(datafolder, "/TargetCells_Bioclim+soils+matches_WY_interpolation.csv"))
