%Script that says how many "neurons" are necessary to fill a voxel

%We test many different numbers of neurons, and for each determine the
%density of neurons per voxel. By fitting a line to the number of
%neurons/voxel as a function of neurons, we can find the number of neurons
%necessary to achieve a specific neuronal density.

neurs=[100:100:2000 100:100:2000 100:100:2000 100:100:2000 100:100:2000];

x=20;
y=20;
z=20;
num_vox=x*y*z;

neuronpervox=zeros(1,length(neurs)); %Vector of how many neurons per voxel

for j=1:length(neurs)
       
num_neur=neurs(j); %number of "neurons"
    

%The following code finds what neurons (1 to n) are in what grid spaces
%Essentially X is a binary matrix with neuron/gridspace matches = 1
%Random center points of voxels
x0=x*rand(num_vox,1);
y0=y*rand(num_vox,1);
z0=z*rand(num_vox,1);


% Simulate "neurons" in voxels

X=zeros(num_neur,num_vox); %Matrix of "neurons" by grid spaces
%The following code finds what neurons (1 to num_neur) are in what grid spaces
%Essentially X is a binary matrix with neuron/gridspace matches = 1

%Loop through neurons and see which voxels each neuron goes through
for i=1:num_neur
    %x_neur, y_neur, z_neur are centers of neurons
    x_neur=rand*x;
    y_neur=rand*y;
    z_neur=rand*z;
    %Create random direction of the neuron 
    dx=rand;
    dy=rand;
    dz=rand;
    m1=dy/dx;
    m2=dz/dy;
    m3=dz/dx;
    th1=atan(m1);
    th2=atan(m2);
    th3=atan(m1*m2);
    
    %The "neurons" are lines of length "len", and we say that a neuron is
    %in a voxel if that neuron is within "dia" of the voxel center. 
    dia=2;
    len=2*sqrt(10);
    
    %Find whether neurons are within distance "dia" from the voxel centers in
    %any of the 3 dimensions
    %We do this by first assuming the neurons are infinitely long lines,
    %and seeing whether the line goes within distance "dia" of a voxel
    %center. Then, we account for the length of the neuron by checking
    %whether the neuron center is within distance "len"/2 from the voxel
    %center.
    
    %Find the distances of all the voxels to the given neuron in each
    %dimension
    deltaX=x0-x_neur;
    deltaY=y0-y_neur;
    deltaZ=z0-z_neur;
    dist_from_line1=abs((deltaY)-m1*(deltaX));
    dist_from_line2=abs((deltaZ)-m2*(deltaY));
    dist_from_line3=abs((deltaZ)-m3*(deltaX));
    %Find the distance from the voxel center to the neuron center
    dist_from_pt=sqrt((deltaX).^2+(deltaY).^2+(deltaZ).^2);
    %Find all the voxels that are close enough to the given neuron 
    q=(dist_from_line1<dia/cos(th1) & dist_from_line2<dia/cos(th2) & dist_from_line3<dia/cos(th3) & dist_from_pt<len/2);
    
    %If a "neuron" is "in" the voxel, put a 1 in the corresponding matrix
    %location
    X(i,q)=1;
end

neuronpervox(j)=nnz(X)/num_vox; %Neuronal density (how many neurons per voxel)
end

%Fit a line to determine an equation for how many neurons/voxel there are
%as a function of the number of neurons. This will be used to determine how
%many neurons are needed to completely fill a voxel (in neur_per_vox_func.m)
[b,bint,r,rint,stats]=regress(neuronpervox',neurs');

%Answer is slope=.0085



