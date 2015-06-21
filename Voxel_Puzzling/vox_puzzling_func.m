function [ r2, err_mean, removed ] = vox_puzzling_func(x,y,z,vox_size,num_vox_removed,method,landmark_pts,plot_figs)
%This is the main function for voxel puzzling, and does all the steps (see
%individual sections)


num_vox=round(x*y*z); %number of voxels

%% Simulate

%Simulate the matrix, Xs, giving which neurons go through which voxels
[ x0, y0, z0, Xs ] = simulate_neur_and_vox(vox_size,x,y,z,num_vox );
%x0,y0,z0 are the x,y,z coordinates of the centers of the voxels

%% Make Similarity Matrix
% clearvars -except X x y n x0 y0 z0 remove

H=transpose(Xs)*Xs;

%Do thresholding (description of methods in paper)

%Thresholding before ULI
if method==1
    thresh=floor(max(H)/2);
    threshmat=repmat(thresh,length(H),1);
    HnnS=H>threshmat;
    %     HnnS=HnnS'; %to go with dist_from_pt2
    HnnS=HnnS+HnnS';
end

%Thresholding before SDM
if method==2
    thresh=min(max(H))-1;
    if thresh==-1
        thresh=0;
    end
    HnnS=H>thresh;
    W=HnnS;
end


%% Remove Voxels

%Remove voxels that aren't part of the main connected component,
%and also remove voxels for testing purposes

vox_remaining=1:num_vox; %Keep track of voxels that remain

%Remove voxels for testing purposes
if num_vox_removed~=0
    shuffle=randperm(num_vox);
    remove_idx=shuffle(1:num_vox_removed);
    HnnS(remove_idx,:)=[]; %Do removal from similarity matrix
    HnnS(:,remove_idx)=[];
    vox_remaining(remove_idx)=[]; %Keep track of voxels that remain
end

%Remove voxels that are not part of the main connected component
[sci, sizes]=scomponents(HnnS); %Find voxels that are part of main connected component
[~,maxidx]=max(sizes);
remove_idx2=find(sci~=maxidx);
HnnS(remove_idx2,:)=[];
HnnS(:,remove_idx2)=[];
vox_remaining(remove_idx2)=[];

removed=num_vox-length(vox_remaining); %Keep track of how many have been removed

% Update the voxel coordinates to not include those that were
% removed
x0_removed=x0(vox_remaining);
y0_removed=y0(vox_remaining);
z0_removed=z0(vox_remaining);

%% Dimensionality reduction

%Use dimensionality reduction techniques on the similarity matrix
%to reconstruct the voxel positions

if method==1 %ULI
    [rec_pos,pts]=uli_func(HnnS,landmark_pts,3);    %Output is the reconstructed position, and the positions of the landmark points (pts)
end

if method==2 %SDM
    [ rec_pos ] = sdm_func( HnnS ); %Output is the reconstructed position
end

%% Test Reconstruction Accuracy (Using R-value and distance error)

%Get correlations between pairwise distances
true_pos=[x0_removed(:) y0_removed(:) z0_removed(:)]; %Matrix of the true positions of the voxel centers
true_dists=pdist(true_pos); %Pairwise distances between all the voxel centers
model_dists=pdist(rec_pos); %Pairwise distances between the reconstructed voxel centers
r2=corr(true_dists',model_dists'); %Correlation between the true distances and the reconstructed distances; note that this is the r-value (not r^2), even though the variable name is r2

%Get reconstruction error (in terms of mean distance)
if ~isnan(r2) %In case there's some problem with the reconstruction (happened in previous versions, but should not anymore)
    
    %Transform the true positions, so the center voxel is at (0,0,0)
    x0_trans=x0_removed-mean(x0_removed);
    y0_trans=y0_removed-mean(y0_removed);
    z0_trans=z0_removed-mean(z0_removed);
    true_pos_trans=[x0_trans y0_trans z0_trans];
    
    %Scale the reconstructed voxel positions so they're on the same
    %scale as the true distances
    scale=mean(true_dists)/mean(model_dists);
    rec_pos_trans=rec_pos*scale; %scaled reconstructed voxel positions
    
    %Center the reconstructed positions
    mean_rec=mean(rec_pos_trans);
    rec_pos_trans(:,1)=rec_pos_trans(:,1)-mean_rec(1);
    rec_pos_trans(:,2)=rec_pos_trans(:,2)-mean_rec(2);
    rec_pos_trans(:,3)=rec_pos_trans(:,3)-mean_rec(3);
    
    %Rotate reconstructed voxels to match with true voxel positions
    [ rec_pos_trans ] = rotate_match_func( rec_pos_trans, true_pos_trans );
    
    %Calculate the distance error
    errmat=[x0_trans-rec_pos_trans(:,1), y0_trans-rec_pos_trans(:,2), z0_trans-rec_pos_trans(:,3)];
    err=sqrt(errmat(:,1).^2+errmat(:,2).^2+errmat(:,3).^2);
    
    err_mean=mean(err); %Save the mean distance error for this loop
    
end

%% Plot Voxels (Original and Reconstructed)

if plot_figs
    
    %Make color map
    xVal=x0/x;
    yVal=y0/x;
    xVal=xVal(:);
    yVal=yVal(:);
    colors_xy = [.6*ones(size(xVal)),xVal,yVal];
    colors_xy=colors_xy(vox_remaining,:);
    
    %Plot the true voxel positions
    figure; scatter3(true_pos_trans(:,1),true_pos_trans(:,2),true_pos_trans(:,3),[],colors_xy,'fill','LineWidth',1)
    title('True Voxel Positions')
    
    %Plot the reconstructed voxel positions
    figure; scatter3(rec_pos_trans(:,1),rec_pos_trans(:,2),rec_pos_trans(:,3),[],colors_xy,'fill','LineWidth',1)
    title('Reconstructed Voxel Positions')
    
    
    % Plot a slice through the reconstructed voxel positions
    for i=10:10
        idx=find(y0>(i-1) & y0<i);
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
    text(.8*max(model_dists),.2*max(true_dists),['r=' num2str(r2)],'FontWeight','bold','FontSize',10)
    box off;
    xlabel('Reconstructed Distances (AU)')
    ylabel('True Distances (voxels)')
    
end