# Puzzle_Imaging
Code for Puzzle Imaging Manuscript

A very brief description of the files is below. Further comments are found within the files.
If you have questions, please email j-glaser@u.northwestern.edu

## Voxel_Puzzling

* **vox_puzzling_main.m**
	* This script does all of the voxel puzzling demonstrations from the manuscript (plotting and accuracy w/ different parameters)
	* Uses vox_puzzling_func.m

* **vox_puzzling_func.m**
	* This is the main function for voxel puzzling, and does all the steps
	* Uses simulate_neur_and_vox.m, uli_func.m, sdm_func.m, rotate_match_func.m, scomponents.m

* **simulate_neur_and_vox.m**
	* This function runs the basic simulation to determine where the voxels are, and which neurons are in which voxels.
	* Uses neur_per_vox_func.m

* **neur_per_vox_func.m**
	* For a vector of voxel sizes, says how many neurons there should be to completely fill the voxels. 
	* Based off results of neurpervox.m

* **neurpervox.m**
	* Script that says how many neurons are necessary to fill a voxel

* **plot_fig3.m**
	* Script to make results plots in Figure 3

* **Data_files**
	* Contains the results used for plotting


## Connectomics_Puzzling

* **conn_puzzling_main.m**
	* This script will do all of the connectomics puzzling applications from the manuscript.
	* Uses conn_puzzling_func.m

* **conn_puzzling_func.m**
	* This is the main connectomics puzzling function
	* Uses simulate_neur_and_conn, simulate_neur_and_conn_layers, sdm_func.m, rotate_match_func.m, scomponents.m

* **simulate_neur_and_conn.m**
	* This function simulates neuron locations and the connections between neurons given a single connection probability function (a single layer)

* **simulate_neur_and_conn_layers.m**
	* This function simulates neuron locations and the connections between neurons given multiple connection probability functions (two layers)

* **plot_fig5.m**
	* Script to make results plots in Figure 3

* **Data_files**
	* Contains the results used for plotting


## Chemical_Puzzling

* **chem_puzzling_main.m**
	* This script produces all the figure outputs for the chemical puzzling applications
	* Uses chem_puzzling_func.m
	* Uses the simulations produced from chem_sim.m

* **chem_puzzling_func.m**
	* This is the main function for chemical puzzling.
	* Uses uli_func.m, rotate_func_2d, scomponents.m

* **chem_sim.m**
	* This function does the chemical puzzling simulation that the results in the paper come from.
	* Uses P3.jpg, the image representing the chemical concentration of the plate

* **rotate_func_2d.m**
	* This function rotates the reconstructed position and gives the sum squared error between this rotation and the true position. It works for 2 dimensional inputs.

* **Data_files**
	* Contains the simulation outputs


## General_Files

* **uli_func.m**
	* This function does the ULI (unweighted landmark isomap) algorithm to reduce dimensions
	* Uses dist_from_pt2.m

* **dist_from_pt2.m**
	* This function outputs the geodesic distance of all points to a given landmark point

* **sdm_func.m**
	* This function does the SDM (sparse diffusion maps) algorithm to reduce dimensions

* **rotate_match_func.m**
	* This function finds the rotation (and potential reflection) to best match the reconstructed position to the true position. 
	* Uses rotate_func.m

* **rotate_func.m**
	* This function rotates the reconstructed position and gives the sum squared error between this rotation and the true position.

* **scomponents.m**
	* Finds the connected components of a graph
	* Uses sparse_to_csr.m
