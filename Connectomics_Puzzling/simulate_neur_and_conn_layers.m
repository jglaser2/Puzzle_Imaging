function [ C, neur_x_loc, neur_y_loc, neur_z_loc ] = simulate_neur_and_conn_layers( x,y,z,scale,p33,p22,p23 )
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

n=x*y*z;

neur_x=x*rand(n,1);
neur_x_loc=scale*neur_x;
neur_y=y*rand(n,1);
neur_y_loc=scale*neur_y;
neur_z=z*rand(n,1);
neur_z_loc=scale*neur_z;

true_pos=[neur_x_loc(:) neur_y_loc(:) neur_z_loc(:)];
true_dist_vec=pdist(true_pos);

layer2=neur_x<x/2;
layer3=neur_x>=x/2;
layers=2*layer2+3*layer3;
layers_dist_vec=pdist(layers); %1 for 2-3 and 0 for 2-2 and 3-3
layers_idx=find(layers_dist_vec);

l2_idx=find(layer2);
l2=zeros(size(layer2));
l2(l2_idx)=rand(size(l2_idx));
l2_dist_vec=pdist(l2);
layer2_idx=find(l2_dist_vec); %Nonzero for 2-2 and 2-3, and zero for 3-3

l3_idx=find(layer3);
l3=zeros(size(layer3));
l3(l3_idx)=rand(size(l3_idx));
l3_dist_vec=pdist(l3);
layer3_idx=find(l3_dist_vec); %Nonzero for 3-3 and 2-3, and zero for 2-2



C_vec=zeros(size(true_dist_vec));
C_vec(layer2_idx)=rand(1,length(layer2_idx))<p22(true_dist_vec(layer2_idx));
C_vec(layer3_idx)=rand(1,length(layer3_idx))<p33(true_dist_vec(layer3_idx));
C_vec(layers_idx)=rand(1,length(layers_idx))<p23(true_dist_vec(layers_idx));

C=squareform(C_vec);

C=C+eye(n);



end

