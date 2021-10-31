#' Find solution to k-points algorithm
#'
#' Determine a number (k) of points that maximize the areal coverage of a study
#' area using a set of matching variables to determine similarity among sites.
#'
#' @param matchingvars Data frame generated using \code{\link{makeInputdata}} or formatted
#' such that: rownames are 'cellnumbers' extracted using the \code{\link{raster::extract}}
#' function, columns 2 and 3 correspond to x and y coordinates, and additional
#' columns correspond to potential matching variables extracted using the
#' \code{\link{raster::rasterToPoints}} function.
#'
#' @param criteria Single value or vector of length equal to the number of matching variables,
#' where values corresponds to the matching criterion for each matching variable
#' in x. If a single value, this will be used as matching criteria for all variables.
#' Default value is 1, corresponding to using raw data for kpoints algorithm.
#'
#' @param klist Single value or vector of values of k to find solutions for.
#' Default value is 200.
#'
#' @param mindelta.area Minimum value for change in area represented between
#' iterations. If the change in area represented is at or below mindelta.area
#' for five consecutive iterations, the kpoints algorithm will assume convergence
#' on an optimal solution and stop. Default value is 50.
#'
#' @param iter Maximum number of iterations before the kpoints algorithm will
#' stop. This parameter prevents kpoints from searching indefinitely for a
#' solution. Default value is 50.
#'
#' @param n_starts The number of random starts (k randomly selected points).
#' For determining the optimal number of points, a small value (e.g., 10) should
#' be sufficient. For finding the final solution for the desired number of points,
#' a larger number (e.g., 100) should be used. Default value is 10.
#'
#' @param areas One of the raster layers used for input data.
#' See \code{\link{raster::area}}.
#'
#' @param verify.stop Boolean. Indicates whether algorithm should display figures
#' to evaluate stopping criteria. Displays a plot of areal coverage vs iteration
#' for each of n_starts. Default is FALSE.
#'
#' @return
#'
#' @export
#' @examples
#' # Load bioclim data for Target Cells (from rMultivariateMatchingAlgorithms package)
#' data(bioclim)
#' # Create data frame of potential matching variables for Target Cells
#' y <- makeInputdata(bioclim)


data(bioclim)
# Create data frame of potential matching variables for Target Cells
allvars <- makeInputdata(bioclim)
#' # Restrict data to matching variables of interest
matchingvars <- allvars[,c("cellnumbers","x","y","bioclim_01","bioclim_04","bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
# Create vector of matching criteria
criteria <- c(0.7,42,3.3,66,5.4,18.4)
# Designate template for areas raster
areas = bioclim[[1]]

# Verify stopping criteria for 200 points
results1 <- kpoints(matchingvars,criteria = criteria,klist = 200,n_starts = 10,mindelta.area = 50,iter = 50,areas = bioclim[[1]],verify.stop = TRUE)


kpoints <- function(matchingvars,criteria = 1,klist = 200,n_starts = 10,mindelta.area = 50,iter = 50,areas,verify.stop = FALSE){
  if (length(criteria == 1)){
    criteria = rep(criteria, (ncol(matchingvars)-3))
  }
  # Calculate areas from template raster
  areas <- raster::area(bioclim[[1]])
  # Standardize variables of interest
  stdvars <- matchingvars[,c("x","y")]
  for (i in 4:ncol(matchingvars)){
       stdvars <- cbind(stdvars, matchingvars[,i]/criteria[i-3])
  }
  # fix colnames
  colnames(stdvars) <- c("x","y",colnames(matchingvars)[4:ncol(matchingvars)])

  # Find total area of study area by extracting areas using cellnumbers
  totalarea <- round(sum(extract(areas, as.numeric(row.names(stdvars)))))

  # Calculate distance matrix for all cells:
  xdist <- distances(stdvars[,3:ncol(stdvars)], id_variable = row.names(stdvars))
  cell_numbers <- rownames(stdvars)


  # Create empty list store k-points solutions and vector to store areas for each k:
  best_list <- list()
  best_area <- vector()

  # run k-points algorithm
  for (k in klist){
    area_history.k <- vector() # store areas from each iteration for each k
    start_list <- list() # store best solution from each start
    # Use area_iterations to verify stopping criteria
    if (verify.stop == TRUE){
    area_iterations <- list() # list of area/iteration/n.start
    }
    # Starts:
    for (st in 1:n_starts){
      # Create vectors to hold history
      point_history <- vector(iter, mode="list")
      area_history <- vector()
      # Select k cells to start, by randomly sampling k cells:
      centers <- stdvars[sample(length(stdvars[,1]), size = k),]
      # Run through i iterations:
      for(i in 1:iter) { # Number of iterations to readjust subset cells:
      print(paste0(Sys.time()," Now starting iteration ",i, " for k = ", k, " of start ", st))
      # Find neighbors
      neighbors <- nearest_neighbor_search(xdist, k = 1, search_indices = which(rownames(stdvars) %in% rownames(centers)))
      # Calculate area within min.dist:
      neigh <- cell_numbers[t(neighbors)]
      stdvars2 <- cbind(stdvars,neigh)
      coverage <- as.numeric(abs(stdvars2[as.character(stdvars2[,ncol(stdvars2)]),3]-stdvars[,3]) <= 1)
      for (cv in 2:6){
        coverage <- cbind(coverage,
                          as.numeric(abs(stdvars2[as.character(stdvars2[,ncol(stdvars2)]),2+cv]-stdvars[,2+cv]) <= 1))
        }
      sum6 <- apply(coverage,1, sum)
      area_history[i] <- round(sum(extract(areas, as.numeric(rownames(stdvars[sum6 == 6,])))))
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
      # Compute cells to represent center of each group:
      # Find centroids of the groups assigned to each of the k-points
      mean_centers <- apply(stdvars, 2, tapply, group_hist[,1], mean)
      # If k > 1, find the actual cells closest to centroids:
      if (k > 1){
        # Add centroids onto bioclim:
        rownames(mean_centers) <- c((nrow(stdvars)+1):(nrow(stdvars)+k))
        stdvarsx <- rbind(stdvars, mean_centers)
        # Compute distance matrix and find nn's to each of centroids.
        xdist1 <- distances(stdvarsx[,3:ncol(stdvarsx)], id_variable = row.names(stdvarsx))
        # Search_indices should be all rows--get nearest two neighbors in case of repeats
        mean_neighbors <- nearest_neighbor_search(xdist1, k = 2, search_indices = c(1:nrow(stdvars)),
                                                  query_indices = c((nrow(stdvars)+1):(nrow(stdvars)+k)))
        mean_neighbors1 <- unique(mean_neighbors[1,])
        # Take care of duplicates (where one point is closest to centroid of two groups)
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
        centers <- stdvars[mean_neighbors1,]
      } else { break }
    }
    # Which iteration had max(area)
    best_run <- which.max(area_history)
    # Insert solution for best iteration for this k, this start into list:
    start_list[[st]] <-  list(center_points=point_history[[best_run]][,1:2], area = area_history[best_run])
    # Record area of best run for this start:
    area_history.k[st] <- area_history[best_run]
    if (verify.stop == TRUE){
    # Add last area_history onto list:
    area_iterations[[st]] <- area_history
    # Plot area_iterations to check stopping criteria
    # Need to verify that it is usually mindelta.area that stops iterations, not iter
    par(mar = c(3,3,2,1), mgp = c(1.5,0.3,0), tcl = -0.2)
    plot(seq(0.5,1, by = 0.1)~c(seq(0,50, length.out = 6)), col = "white", xlab = "Number of Iterations",
         ylab = "Proportion of area covered", ylim = c(0,1), cex.lab = 1.2, cex.axis = 1.2,
         main = paste0("Verify stopping criteria: k = ", k))
    for (i in 1:length(area_iterations)){
      lines(c(area_iterations[[i]]/totalarea)~c(1:length(area_iterations[[i]])), col = rainbow(length(area_iterations))[i], type = "l",
            lwd = 1.5)
    }
    }
    }
  # Insert solution and area for each k into best_list and best_area lists
  best_list[[which(klist == k)]] <- start_list[[which.max(area_history.k)]]$center_points
  names(best_list)[which(klist==k)] <- paste0("k_",k)
  best_area[[which(klist == k)]] <- start_list[[which.max(area_history.k)]]$area
  names(best_area)[which(klist==k)] <- paste0("k_",k)
}
# Create list of output results
results <- list()
results$solutions <- best_list
results$solution_areas <- best_area
results$totalarea <- totalarea
results$klist <- klist
results$iter <- iter
results$n_starts <- n_starts
results$criteria <- criteria
results$mindelta.area <- mindelta.area
return(results)
}


#' Determine optimal number of points (k)
#'
#' Determine a number (k) of points that maximize the areal coverage of a study
#' area using a set of matching variables to determine similarity among sites.
#'
#' @param x Output from \code{\link{kpoints}} function.
#'
#' @return Plot of the proportion of the study area covered for each value of k,
#' or if only one value of k was used, reports coverage for that solution.
#'
#' @export
#' @examples
#' # Load bioclim data for Target Cells (from rMultivariateMatchingAlgorithms package)
#' data(bioclim)
#' Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(bioclim)
#' # Subset to include only matching variables
#' matchingvars <- allvars[,c("cellnumbers","x","y","bioclim_01","bioclim_04","bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
#' # Create vector of matching criteria
#' criteria <- c(0.7,42,3.3,66,5.4,18.4)
#' # Designate template for areas raster
#' areas = bioclim[[1]]
#' # Create sequence of values for k
#' klist = seq(25,100, by = 25)
#' # Run kpoints algorithm for klist
#' results3 <- kpoints(matchingvars,criteria = criteria,klist = klist,n_starts = 1,mindelta.area = 50,iter = 50,areas = bioclim[[1]])
#' # Find optimal number of points (k)
#' plotcoverage(results3)

plotcoverage <- function(x){
  if (length(x["solution_areas"][[1]])>1){
par(tcl = 0.3, mgp = c(1.5,0.2,0), mar = c(3,3,1,1))
plot(x["solution_areas"][[1]]/x["totalarea"][[1]]~x["klist"][[1]],
     ylim = c(0,1), type = "b",
     main = "Coverage for different numbers of k-points",
     ylab = "Proportion of area represented",
     xlab = expression(italic("Number of points (k)")),
     cex.lab = 1.2, cex.axis = 1.2, pch = 16, lwd = 2)
  } else {
    print(paste0("Coverage for k = ",x["klist"][[1]]," is ", x["solution_areas"][[1]], " km2, (",
                 round(x["solution_areas"][[1]]/x["totalarea"][[1]]*100),"%)."))
  }
}
