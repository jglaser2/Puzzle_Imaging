function [ r2_vec, err_mean, removed ] = conn_puzzling_func( x,y,z,scale,plot_flag,p33,p22,p23 )
%This is the main connectomics puzzling function

%This function takes as input:
%-the number of neurons that fit going in the x,y,and z directions (x, y, and z)
%-how far apart the neurons are from each other in um (scale)
%-probability distributions of connections as a function of distance between neurons (p33, p22, p23 - see next comment)
%-whether to plot (plot_flag)

%If there is one layer, there is one probability distribution (p33)
%If there are two layers, then there can be three probability
%distributions: p33 and p22 (probability functions w/in the two layers) and
%p23 (probability function between the two layers).

%As output, it gives:
%-the correlation between the true and reconstructed pairwise distances between neurons (r2_vec)
%-the average distance error between the true and reconstructed neuron locations (err_mean)
%-the number of neurons that needed to be excluded because they were not connected to the main connected component of the similarity matrix (removed).

narginchk(6,8);

%% Simulate Neurons and Connections Randomly

%Simulate the matrix,C, giving which neurons are connected to each other
if nargin==6
    [ C, neur_x_loc, neur_y_loc, neur_z_loc ]=simulate_neur_and_conn( x,y,z,scale,p33 );
end
if nargin==8
    [ C, neur_x_loc, neur_y_loc, neur_z_loc ]=simulate_neur_and_conn_layers( x,y,z,scale,p33,p22,p23 );
end
%neur_x_loc, neur_y_loc, and neur_z_loc are the x,y,z coordinates of the centers of the voxels

%% Remove neurons that are not part of the main connected component

%Note that this section is most likely not necessary.
%In my simulations, there were generally no neurons that had to be removed

[sci, sizes]=scomponents(C); %Find neurons that are part of main connected component
[~,maxidx]=max(sizes);
remove_idx2=find(sci~=maxidx);
C(remove_idx2,:)=[]; %Do removal from similarity (connectivity) matrix
C(:,remove_idx2)=[];

removed=length(remove_idx2); %Keep track of how many have been removed

% Update the neuron coordinates to not include those that were
% removed
neur_x_loc_removed=neur_x_loc;
neur_x_loc_removed(remove_idx2)=[];
neur_y_loc_removed=neur_y_loc;
neur_y_loc_removed(remove_idx2)=[];
neur_z_loc_removed=neur_z_loc;
neur_z_loc_removed(remove_idx2)=[];

%% Dimensionality Reduction

%Use dimensionality reduction techniques on the similarity matrix
%to reconstruct the neuron positions
[ rec_pos ] = sdm_func( C );


%% Test Reconstruction Accuracy (Using R-value and distance error)

%Get correlations between pairwise distances

true_pos=[neur_x_loc_removed(:) neur_y_loc_removed(:) neur_z_loc_removed(:)]; %Matrix of the true positions of the neurons
true_dists=pdist(true_pos); %Pairwise distances between all the neurons
model_dists=pdist(rec_pos); %Pairwise distances between the reconstructed neuron locations
r2_vec=corr(true_dists',model_dists');

%Get reconstruction error (in terms of mean distance)

%Transform the true positions, so the center neuron is at (0,0,0)
neur_x_trans=neur_x_loc_removed-mean(neur_x_loc_removed);
neur_y_trans=neur_y_loc_removed-mean(neur_y_loc_removed);
neur_z_trans=neur_z_loc_removed-mean(neur_z_loc_removed);
neur_loc=[neur_x_trans neur_y_trans neur_z_trans];

%Scale the reconstructed neuron positions so they're on the same
%scale as the true distances
scale2=mean(true_dists)/mean(model_dists);
rec_pos_trans=rec_pos*scale2;

%Center the reconstructed positions
mean_rec=mean(rec_pos_trans);
rec_pos_trans(:,1)=rec_pos_trans(:,1)-mean_rec(1);
rec_pos_trans(:,2)=rec_pos_trans(:,2)-mean_rec(2);
rec_pos_trans(:,3)=rec_pos_trans(:,3)-mean_rec(3);

%Rotate reconstructed neurons to match with true neuron positions
[ rec_pos_trans ] = rotate_match_func( rec_pos_trans, neur_loc );

%Calculate the distance error
errmat=[neur_x_trans-rec_pos_trans(:,1), neur_y_trans-rec_pos_trans(:,2), neur_z_trans-rec_pos_trans(:,3)];
err=sqrt(errmat(:,1).^2+errmat(:,2).^2+errmat(:,3).^2); %Distance error of every neuron
err_mean=mean(err);


%% Plot Reconstruction

if plot_flag
    
    %Make color map
    xVal=neur_x_loc_removed/x/scale;
    yVal=neur_y_loc_removed/y/scale;
    colors_xy = [.6*ones(size(xVal)),xVal,yVal];
    %Plot the reconstructed neuron positions
    figure; scatter3(rec_pos_trans(:,1),rec_pos_trans(:,2),rec_pos_trans(:,3),[],colors_xy,'fill','LineWidth',1)
    title('Reconstructed Voxel Positions')
    
    
    % Plot a slice through the reconstructed neuron positions
    for i=10:10
        idx=find(neur_y_loc>scale*(i-1) & neur_y_loc<scale*i);
        figure; scatter3(rec_pos_trans(idx,1),rec_pos_trans(idx,2),rec_pos_trans(idx,3),[],colors_xy(idx,:),'fill','LineWidth',1)
    end
    set(gca,'Visible','off')
    title('Reconstructed Voxel Positions Slice')
    
    
    
    
    %Plot pairwise model (reconstructed) distances against pairwise true distances
    %Essentially, what the below code does is bins all the
    %modelpairwise distances, and then finds the mean true pairwise
    %distances for all voxel pairs in that bin
    
    %Make bins based on all model pairwise distances
    num_bins=100;
    max_dist=max(model_dists);
    bin_dist=max_dist/num_bins;
    bin_vec=(1:num_bins)*bin_dist;
    
    %Find true pairwise distances in each of those bins
    m=zeros(1,num_bins); %mean pairwise dists
    s=zeros(1,num_bins); %std pairwise dists
    for i=1:num_bins
        idx=find(model_dists<bin_dist*i & model_dists>bin_dist*(i-1));
        m(i)=mean(true_dists(idx));
        s(i)=std(true_dists(idx));
    end
    
    figure; %Plot
    plot(bin_vec,m);
    hold on;
    plot(bin_vec,m-s,'g');
    plot(bin_vec,m+s,'g');
    text(.8*max(model_dists),.2*max(true_dists),['r=' num2str(r2_vec)],'FontWeight','bold','FontSize',10)
    box off;
    xlabel('Reconstructed Distances (AU)')
    ylabel('True Distances (um)')
    
    
    
end

end
%Note: the below were used to see how reconstruction was within individual
%layers. To use this, layer2 and layer3 need to be made outputs of
%"simulate_neur_and_conn_layers.m"

% model_dists2=pdist(rec_pos(layer2,1:dim));
% true_dists2=pdist(true_pos(layer2,:));
% r2_2=corr(true_dists2',model_dists2');
%
% model_dists3=pdist(rec_pos(layer3,1:dim));
% true_dists3=pdist(true_pos(layer3,:));
% r2_3=corr(true_dists3',model_dists3');
