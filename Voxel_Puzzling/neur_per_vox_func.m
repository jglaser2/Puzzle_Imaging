function [ neurs ] = neur_per_vox_func( vox_size )
%For a vector of voxel sizes (vox_size), says how many neurons (neurs) there should be to
%completely fill the voxels. 

slope = .0085; 
%This is the slope of the graph of the number of neurons/voxel as a function of neurons
%which was calculated in "neurpervox.m"

%We want to calculate the number of neurons of cross-sectional area 1 um^2
%that we would need so that each voxel with a side length of vox_size^2 um^2 
%would be filled, assuming the neurons go completely through the voxels (see manuscript for more details)

%We need vox_size^2 neurons / voxel
neurs=round(vox_size.^2/slope);

end

