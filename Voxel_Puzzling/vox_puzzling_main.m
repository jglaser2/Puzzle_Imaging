%This script does all of the voxel puzzling demonstrations from the manuscript (plotting and accuracy w/ different parameters)
%Just run the cell for the purpose you want

%% Plots (Initial and reconstructed voxels and r2 plot - Fig 3A&B)

method=1; %Dimensionality Reduction method: ULI=1, %SDM=2
landmark_pts=10; %Number of landmark points to use in ULI method.
plot_figs=1;

vox_size=5; %Size of voxels (average length of a side, in um)
num_vox_removed=0;

%SIMULATION SIZE
x=20; %Size of simulation (x dimension; in number of voxel average edge lengths)
y=20; %Size of simulation (y dimension; in number of voxel average edge lengths)
z=20; %Size of simulation (z dimension; in number of voxel average edge lengths)

vox_puzzling_func(x,y,z,vox_size,num_vox_removed,method,landmark_pts,plot_figs);


%% Test voxel size (data for Fig 3C,D,E)

num=200; %number of simulations -- (used for test==1 and test==2)
vox_size_vec=2:10; %Vector of size of voxels to test (average side lengths in um)
num_vox_removed=0;

method=1; %Dimensionality Reduction method: ULI=1, %SDM=2
landmark_pts=10; %Number of landmark points to use in ULI method.

plot_figs=0;

x=20; %Size of simulation (x dimension; in number of voxel average edge lengths)
y=20; %Size of simulation (y dimension; in number of voxel average edge lengths)
z=20; %Size of simulation (z dimension; in number of voxel average edge lengths)

num_loops=length(vox_size_vec);

r2=zeros(num_loops,num);
err_mean=zeros(num_loops,num);
removed=zeros(num_loops,num);

for jj=1:num_loops %Loop over different voxel sizes or numbers of voxels removed
    
    vox_size=vox_size_vec(jj);
    for kk=1:num %Loop over different simulations/reconstructions
        
        
        [ r2(jj,kk), err_mean(jj,kk), removed(jj,kk) ] = vox_puzzling_func(x,y,z,vox_size,num_vox_removed,method,landmark_pts,plot_figs);
        
        
    end
end


%% Test number of voxels removed (data for Fig 3F,G)

num=200; %number of simulations -- (used for test==1 and test==2)
vox_size=5; 
perc_vox_removed_vec=0:.05:.9; %Vector of percentages of voxels removed to test -- (used for test==2

method=2; %Dimensionality Reduction method: ULI=1, %SDM=2
landmark_pts=10; %Number of landmark points to use in ULI method.

plot_figs=0;

x=20; %Size of simulation (x dimension; in number of voxel average edge lengths)
y=20; %Size of simulation (y dimension; in number of voxel average edge lengths)
z=20; %Size of simulation (z dimension; in number of voxel average edge lengths)
num_vox=x*y*z;

num_loops=length(perc_vox_removed_vec);

r2=zeros(num_loops,num);
err_mean=zeros(num_loops,num);
removed=zeros(num_loops,num);

for jj=1:num_loops %Loop over different voxel sizes or numbers of voxels removed
    
    num_vox_removed=num_vox*perc_vox_removed_vec(jj);
    
    for kk=1:num %Loop over different simulations/reconstructions
                
        [ r2(jj,kk), err_mean(jj,kk), removed(jj,kk) ] = vox_puzzling_func(x,y,z,vox_size,num_vox_removed,method,landmark_pts,plot_figs);
                
    end
end


%% Plot reconstruction of rectangular shape- Fig S2)

method=1; %Dimensionality Reduction method: ULI=1, %SDM=2
landmark_pts=10; %Number of landmark points to use in ULI method.
plot_figs=1;

vox_size=5; %Size of voxels (average length of a side, in um)
num_vox_removed=0;

%SIMULATION SIZE
%Make z half as large as x and y, with x*y*z=8000 (same as with cube)
x=2*2000^(1/3); %Size of simulation (x dimension; in number of voxel average edge lengths)
y=2*2000^(1/3); %Size of simulation (y dimension; in number of voxel average edge lengths)
z=2000^(1/3); %Size of simulation (z dimension; in number of voxel average edge lengths)

vox_puzzling_func(x,y,z,vox_size,num_vox_removed,method,landmark_pts,plot_figs);