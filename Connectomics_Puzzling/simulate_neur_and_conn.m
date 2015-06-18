function [ C, neur_x_loc, neur_y_loc, neur_z_loc ] = simulate_neur_and_conn( x,y,z,scale,p )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

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
% true_dist_matrix=squareform(true_dist_vec);

C_vec=rand(1,length(true_dist_vec))<p(true_dist_vec);
C=squareform(C_vec);
C=C+eye(n);


end

