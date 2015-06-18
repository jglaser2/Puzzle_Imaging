%This script produces all the figure outputs for the chemical puzzling
%applications

%% Plot Reconstruction w/ Cell Density=1%, and 50 Basepairs

load('sim_dens=1_conj=100.mat') %Load the simulation giving the pioneer cell locations and connections
num_landmarks=5; %Number of landmark points for ULI algorithm
bc_length=50; %Number of basepairs used to encode the chemical concentration
plot_toggle=1; %Plot

chem_puzzling_func(dists,bc_length,num_landmarks,plot_toggle)

%% Plot Reconstruction w/ Cell Density=1%, and 2 Basepairs

load('sim_dens=1_conj=100.mat') %Load the simulation giving the pioneer cell locations and connections
num_landmarks=5; %Number of landmark points for ULI algorithm
bc_length=2; %Number of basepairs used to encode the chemical concentration
plot_toggle=1; %Plot

chem_puzzling_func(dists,bc_length,num_landmarks,plot_toggle);

%% Plot Reconstruction w/ Cell Density=0.1%, and 50 Basepairs

load('sim_dens=pt1_conj=100.mat') %Load the simulation giving the pioneer cell locations and connections
num_landmarks=5; %Number of landmark points for ULI algorithm
bc_length=50; %Number of basepairs used to encode the chemical concentration
plot_toggle=1; %Plot

chem_puzzling_func(dists,bc_length,num_landmarks,plot_toggle)

%% Plot Reconstruction w/ Cell Density=0.1%, and 2 Basepairs

load('sim_dens=pt1_conj=100.mat') %Load the simulation giving the pioneer cell locations and connections
num_landmarks=5; %Number of landmark points for ULI algorithm
bc_length=2; %Number of basepairs used to encode the chemical concentration
plot_toggle=1; %Plot

chem_puzzling_func(dists,bc_length,num_landmarks,plot_toggle)

%% Plot Reconstruction w/ Conjugation Probability = 40% (Cell Density=1%, 50 basepairs)

load('sim_dens=1_conj=40.mat') %Load the simulation giving the pioneer cell locations and connections
num_landmarks=5; %Number of landmark points for ULI algorithm
bc_length=50; %Number of basepairs used to encode the chemical concentration
plot_toggle=1; %Plot

chem_puzzling_func(dists,bc_length,num_landmarks,plot_toggle)

%% Plot Reconstruction w/ Conjugation Probability = 30% (Cell Density=1%, 50 basepairs)

load('sim_dens=1_conj=30.mat') %Load the simulation giving the pioneer cell locations and connections
num_landmarks=5; %Number of landmark points for ULI algorithm
bc_length=50; %Number of basepairs used to encode the chemical concentration
plot_toggle=1; %Plot

chem_puzzling_func(dists,bc_length,num_landmarks,plot_toggle)

%% Plot Reconstruction w/ Conjugation Probability = 20% (Cell Density=1%, 50 basepairs)

load('sim_dens=1_conj=20.mat') %Load the simulation giving the pioneer cell locations and connections
num_landmarks=5; %Number of landmark points for ULI algorithm
bc_length=50; %Number of basepairs used to encode the chemical concentration
plot_toggle=1; %Plot

chem_puzzling_func(dists,bc_length,num_landmarks,plot_toggle)
