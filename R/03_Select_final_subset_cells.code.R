##########################################
#
# Purpose: Select final subset cells (200)
#
# Summary: Use k-points algorithm to select 200 points for simulation
#
# Inputs: "drylands_wyoming.tif"
#         "TargetCells_Bioclim+soils_WY.csv.csv"
#
# Outputs: "SubsetCells_siteselction.csv"
#
# Rachel Renne
# March 24, 2021
#
##########################################

library(raster)
library(distances)

########################################
# Step 1: Read in data and set up for k-points algorithm

# Folder where input and output data are stored
datafolder <- getwd()

# Read in template raster from which we can generate areas
# for use in determining areal coverage of different sets of subset cells
wydrylands <- raster(paste0(datafolder,"/drylands_wyoming.tif"))
areas <- area(wydrylands) # Gives approx. area in km2
rm(wydrylands)

# Read in bioclim and soils data
bioclim <- readRDS(paste0(datafolder,"/TargetCells_Bioclim+soils_WY.csv"))
# Note: cellnumbers in this dataframe correspond to the unique cellnumbers associated with each
# target cell. These serve as unique identifiers throughout the code and can
# be extracted from the raster datasets that define the target cells using the 
# extract function in the raster package with "cellnumbers" = TRUE

# Find total area of trimmed extent by extracting areas using cellnumbers
totalarea <- round(sum(extract(areas, as.numeric(row.names(bioclim)))))
# 246771 km2

# limit to just the 6 variables of interest (i.e., 6 bioclim matching variables):
bioclim1a <- data.frame(bioclim[,c(2,3,4,7,12,15,18,21)])

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

# Standardize variables of interest according to percentage of range
bioclim1 <- cbind(bioclim1a[,1:2],bioclim1a$bioclim_01/0.7, 
                  bioclim1a$bioclim_04/42, 
                  bioclim1a$bioclim_09/3.3,
                  bioclim1a$bioclim_12/66,
                  bioclim1a$bioclim_15/5.4, 
                  bioclim1a$bioclim_18/18.4)
# fix column names
colnames(bioclim1) <- c(colnames(bioclim1a)[1:2],"bioclim_01","bioclim_04","bioclim_09", "bioclim_12","bioclim_15","bioclim_18")

#########################################################################
# Step 2: Select final subset cells (200) by running k-points algorithm with 100 random starts

# After each solution for k is complete, algorithm will generate a plot of 
# proportion of area covered ~ number of iterations
# This will verify that the areal coverage stopping criterion stops the algorithm: 
# i.e., that the maximum iterations (50 in our example) are not reached and that
# coverage levels out before stopping

# Calculate distance matrix for all cells:
xdist <- distances(bioclim1[,3:ncol(bioclim1)], id_variable = row.names(bioclim1))
cell_numbers <- rownames(bioclim1)

# Settings for kpoints
mindelta.area = 50 # Iterations stop if area gained between iterations is < 50 km for 5 consecutive iterations
iter = 50 # Maximum number of iterations for each random start
n_starts = 100 # Number of random starts (k randomly selected points)
# klist is a single value--selected previously as the optimal value for n:
klist <- 200

# Create empty list store kpoints solution and vector to store areas for each k (just one in this case):
best_list <- list()
best_area <- vector()

# run k-points algorithm
for (k in klist){  
  area_history.k <- vector() # store area for each k
  start_list <- list() # store best solution from each start
  # Use area_iterations to verify stopping criteria
  area_iterations <- list() # list of area/iteration/n.start
  # Starts:
  for (st in 1:n_starts){
    # Create vectors to hold history
    point_history <- vector(iter, mode="list")
    area_history <- vector()
    # Select k cells to start, by randomly sampling k cellss:
    # Can "set seed" for reproducibility, if desired:
    # set.seed(356)
    centers <- bioclim1[sample(length(bioclim1[,1]), size = k),]
    # Run through i iterations:
    for(i in 1:iter) { # Number of iterations to readjust subset cells:
      print(paste0(Sys.time()," Now starting Iteration ",i, " for K = ", k, " of Start ", st))
      # Find neighbors
      neighbors <- nearest_neighbor_search(xdist, k = 1, search_indices = which(rownames(bioclim1) %in% rownames(centers)))
      # Calculate area within min.dist:
      neigh <- cell_numbers[t(neighbors)]
      bioclim2 <- cbind(bioclim1,neigh)
      coverage <- as.numeric(abs(bioclim2[as.character(bioclim2[,ncol(bioclim2)]),3]-bioclim1[,3]) <= 1)
      for (cv in 2:6){
        coverage <- cbind(coverage,
                          as.numeric(abs(bioclim2[as.character(bioclim2[,ncol(bioclim2)]),2+cv]-bioclim1[,2+cv]) <= 1))
      }
      sum6 <- apply(coverage,1, sum)  
      area_history[i] <- round(sum(extract(areas, as.numeric(rownames(bioclim1[sum6 == 6,])))))
      print(paste0("For iteration ",i," with ", k," cells, we cover ", area_history[i]," km2 (",round(area_history[i]/totalarea*100),"%)."))
      group_hist <- cbind(cell_numbers[t(neighbors)], cell_numbers)
      point_history[[i]] <- centers
      # Mindelta.area stopping critera executed below:
      if (i > 5){
        delta.area1 <- area_history[i] - area_history[[i-1]]
        delta.area2 <- area_history[[i-1]]-area_history[[i-2]]
        delta.area3 <- area_history[[i-2]]-area_history[[i-3]]
        delta.area4 <- area_history[[i-3]]-area_history[[i-4]]
        delta.area5 <- area_history[[i-4]]-area_history[[i-5]]
        if (delta.area1 <= mindelta.area &&
            delta.area2 <= mindelta.area &&
            delta.area3 <= mindelta.area &&
            delta.area4 <= mindelta.area &&
            delta.area5 <= mindelta.area){ break }
      }
      # Compute cells to represent center of each group (i.e. group of nn to each of k points):
      # Find centroids of the groups assigned to each of the k-points
      mean_centers <- apply(bioclim1, 2, tapply, group_hist[,1], mean) # Calculate centroid of each group
      # If k > 1, find the actual cells closest to centroids:
      if (k > 1){ 
        # Add centroids onto bioclim:
        rownames(mean_centers) <- c((nrow(bioclim1)+1):(nrow(bioclim1)+k))
        bioclimx1 <- rbind(bioclim1, mean_centers)
        # Compute distance matrix and find nn's to each of centroids.
        xdist1 <- distances(bioclimx1[,3:ncol(bioclimx1)], id_variable = row.names(bioclimx1))
        # Search_indices should be all rows--get nearest two neighbors in case of repeats
        mean_neighbors <- nearest_neighbor_search(xdist1, k = 2, search_indices = c(1:nrow(bioclim1)), 
                                                  query_indices = c((nrow(bioclim1)+1):(nrow(bioclim1)+k)))
        mean_neighbors1 <- unique(mean_neighbors[1,])
        # Take care of duplicates
        if (length(unique(mean_neighbors1)) < k){
          dups <- table(mean_neighbors[1,])
          dups1 <- as.numeric(names(dups[dups > 1]))
          dups2 <- mean_neighbors[,mean_neighbors[1,] %in% dups1]
          prob_len <- ncol(dups2)
          dup3 <- unique(dups2[1,])
          dup4 <- vector()
          for (d in (length(unique(dups2[1,]))+1):prob_len){
            dupx <- dups2[,dups2[1,]==dup3[d-length(unique(dups2[1,]))]]
            if (!is.na(dupx[1,1])){
              for (cc in 2:ncol(dupx)){
                dup4 <- append(dup4, dupx[2,cc])
              } }
          }
          
          mean_neighbors1 <-c(unique(mean_neighbors[1,]),dup4)
        }
        centers <- bioclim1[mean_neighbors1,]
      } else { break }
    }
    # Which iteration had max(area)
    best_run <- which.max(area_history)
    # Insert best iteration for this k, this start into list:
    start_list[[st]] <-  list(center_points=point_history[[best_run]], area = area_history[best_run])
    # Record area of best run for this start:
    area_history.k[st] <- area_history[best_run]
    # Add last area_history onto list:
    area_iterations[[st]] <- area_history
  }
  # Plot area_iterations to check stopping criteria
  # Need to verify that it is usually mindelta.area that stops iterations, not iter
  par(mar = c(3,3,2,1), mgp = c(1.5,0.3,0), tcl = -0.2)
  plot(seq(0.5,1, by = 0.1)~c(seq(0,50, length.out = 6)), col = "white", xlab = "Number of Iterations",
       ylab = "Proportion of area covered", ylim = c(0,1), cex.lab = 1.2, cex.axis = 1.2,
       main = paste0("Verify stopping criteria: k = ", k))
  for (i in 1:length(area_iterations)){
    lines(c(area_iterations[[i]]/totalarea)~c(1:length(area_iterations[[i]])), col = rainbow(length(area_iterations))[i], type = "l",
          lwd = 2)
  }
  best_list[[which(klist == k)]] <- start_list[[which.max(area_history.k)]]
  best_area[[which(klist == k)]] <- start_list[[which.max(area_history.k)]]$area
}

# Values may vary somewhat from run to run due to random selections of cells
# Results from run on 3/24/2021
# best_area <- 223410
# best_area/totalarea = 0.905

# Pull center points out of saved list
center_points <- best_list[[1]]$center_points
# Add a cellnumbers column
center_points$cellnumbers <- row.names(center_points)

# Save file
#write.csv(center_points, paste0(datafolder, "/SubsetCells_siteselction.csv"))
