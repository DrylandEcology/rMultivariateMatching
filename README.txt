Multivariate matching algorithms for site selection (A) and interpolation (B) as described in:
Renne, R.R., Schlaepfer, D.R., Palmquist, K.A., Lauenroth, W.K., & Bradford, J.B. In preparation. 
  Estimating complex ecological variables at high resolution in heterogeneous terrain using a multivariate matching algorithm.

Description:
This repository includes code for interpolating multivariate and time series simulation output to high resolution maps as described by Renne et al. (in prep). 
  First, we apply a variation of k-means clustering to select a set of sites for simulation that maximizes the area represented for a given number 
  of sites (A. Site Selection). Then, we used multivariate matching to interpolate simulation results to high-resolution maps (B. Interpolation). 
  Both methods rely on a user-defined set of matching variables that are assigned weights such that matched sites will be within a prescribed range 
  for each variable. We demonstrate each method with case studies using an individual-based plant simulation model to illustrate site selection (A) 
  and an ecosystem water balance simulation model for interpolation (B). The code presented here is intended to illustrate both methods for drylands 
  in Wyoming, USA. Datasets referenced in the code are available at https://app.box.com/s/rp7zcs633eih6j6a8s9t2hmrkkgudkcf

A. Site Selection 
	I. Select matching variables and matching criteria (code file: "01_Identify_matching_variables_and_criteria.code.R")
		Inputs: "drylands_wyoming.tif"
				"TargetCells_Bioclim+soils_WY.csv"
		Outputs: matching criteria

	II. Determine optimal number of sites to simulate (code file: "02_Determine_number_of_subset_cells.code.R")
		Inputs: "drylands_wyoming.tif"
				"TargetCells_Bioclim+soils_WY.csv"
		Outputs: "Coverage_by_n_siteselection.png"
				 Optimal number of sites (n = 200 in our example)
	
	III. Select final Subset cells with 100 random starts of k-points algorithm (code file: "03_Select_final_subset_cells.code.R")
		Inputs: "drylands_wyoming.tif"
				"TargetCells_Bioclim+soils_WY.csv.csv"
		Outputs: "SubsetCells_siteselction.csv"

	IV. Evaluate matching between Target and Subset cells in three ways for selected sites (code file: "04_Evaluating_matching_siteselection.code.R")
		Inputs: "TargetCells_Bioclim+soils_WY.csv"
				"SubsetCells_siteselction.csv"
				"drylands_wyoming.tif"
		Outputs: "Mapped Coverage_siteselection.png"
				 "Variable Coverage_siteselection.png"
				 "SD_Diffs_allcells+matched_siteselection.png"
				 "Distances_subset+target_cells_siteselection.png"
				 "Distances_adjacent_subset_cells_siteselection.png"
				 
B. Interpolation
	I. Find best match for each Target cells by matching using 2-step process: 1) match to climate variables, 2) match to soils (code file: "05_Two_step_Matching_v3.code.R")
		Inputs: "TargetCells_Bioclim+soils_WY.csv"
				"SubsetCells_bioclim+soils+results_WY_interpolation.csv"
		Outputs: "TargetCells_Bioclim+soils+matches_WY_interpolation.csv"
				 
	II. Evaluate matching between Target and Subset cells in three ways for selected sites (code file: "06_Evaluating_matching_interpolation.code.R")
		Inputs: "drylands_wyoming.tif"
				"TargetCells_Bioclim+soils+matches_WY_interpolation.csv"
		Outputs: "Mapped Coverage_interpolation.png"
				 "Mapped Coverage_soils_interpolation.png"
				 "Variable Coverage_interpolation.png"
				 "SD_Diffs_allcells+matched_interpolation.png"
				 "Distances_subset+target_cells_interpolation.png"
				 "Distances_adjacent_subset_cells_interpolation.png"
				 
	III. Leave-one-out cross validation to estimate matching errors (code file: "07_Leave_one_out_CV_code.R")	
		Inputs: "SubsetCells_bioclim+soils+results_WY_interpolation.csv"
		Outputs: "SD_Diffs_LOOCV_interpolation.png"
				 "Hist_Estimated_Errors_LOOCV_interpolation.png"	
	
	IV. Interpolate simulation results to create high-resolution maps (code file: "08_Interpolating_results_to_highres_maps.code.R")
		Inputs: "SubsetCells_bioclim+soils+results_WY_interpolation.csv"
				"TargetCells_Bioclim+soils+matches_WY_interpolation.csv"
				"drylands_wyoming.tif"
		Outputs: "Interpolated_output_vars_maps.png"