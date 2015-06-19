function [ C, neur_x_loc, neur_y_loc, neur_z_loc ] = simulate_neur_and_conn( x,y,z,scale,p )
%This function simulates neuron locations and the connections between
%neurons given a single connection probability function (a single layer)

%Its inputs are:
%x,y,z - the number of neurons that fit going in the x,y,and z directions
%scale - how far apart the neurons are from each other in um
%p - probability distributions of connections as a function of distance between neurons

%Its outputs are:
%C - the connectivity matrix
%neur_x_loc, neur_y_loc, neur_z_loc - x,y, and z coordinates of the neuron locations


% Simulate Neurons Randomly
n=x*y*z;

neur_x=x*rand(n,1);
neur_x_loc=scale*neur_x;
neur_y=y*rand(n,1);
neur_y_loc=scale*neur_y;
neur_z=z*rand(n,1);
neur_z_loc=scale*neur_z;

% Determine connections

true_pos=[neur_x_loc(:) neur_y_loc(:) neur_z_loc(:)];

true_dist_vec=pdist(true_pos);

C_vec=rand(1,length(true_dist_vec))<p(true_dist_vec);
C=squareform(C_vec);
C=C+eye(n);

end

