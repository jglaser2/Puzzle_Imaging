function [r2] = chem_puzzling_func(dists,bclength,num_landmarks,plottoggle)
%This is the main function for chemical puzzling.

%Inputs:
%dists-the struct that is output from the simulations, which gives the locations of the
%pioneer cells (dists.init), the image/plate being used (dists.images), and
%the connectivity matrix between pioneer cells (dists.josh)
%bclength- Number of basepairs used to encode the chemical concentration
%num_landmarks - number of landmark points used in ULI algorithm
%plottoggle - whether to plot

%Outputs:
%r2- the correlation between pairwise distances between all pioneer cells originally and in reconstruction


init_idxs = [1 2 3];

%% Remove cells that are not part of the main connected component
%Find initial points that need to be removed due to not being part of the
%main connected component
C0=dists.josh;
[sci, sizes]=scomponents(C0); %Find cells that are part of main connected component
[~,maxidx]=max(sizes);
remove=find(sci~=maxidx);

%Remove from the initial locations
init_x=dists.init.x;
init_x(remove)=[];
init_y=dists.init.y;
init_y(remove)=[];

%Remove from the connectivity matrix
C=C0;
C(remove,:)=[];
C(:,remove)=[];
%Also add in self connections
diags=max(C);
C=C+diag(diags);

%% Dimensionality Reduction (Reconstruct Points)

%Get matrix for reconstruction and reconstruct
Cs=sparse(C);
C2=Cs^2>0; % this makes a binary connectivity matrix of two steps -- anything that is two connection steps will be 1, anything that is 3 or more will be 0.
[rec_pos,pts]=uli_func(C2,num_landmarks,2);

%% Rotate the reconstruction (and reflection) to match the true locations based on knowledge of the 3 reference points

%First we need to make sure the initial points have the same scale as the
%true locations, and make sure that everything is centered (for rotation
%matrix purposes)

true_pos_init=[init_x(init_idxs,:) init_y(init_idxs,:)];
true_dists=pdist(true_pos_init); %Pairwise distances between the true locations of the initial pioneer cells
model_dists=pdist(rec_pos(init_idxs,:)); %Pairwise distances between the reconstructed locations of the initial pioneer cells

%Scale the reconstructed cell positions so they're on the same scale as the true distances
scale=mean(true_dists)/mean(model_dists); 
rec_pos_trans_init=rec_pos(init_idxs,:)*scale;

true_pos_trans_init=true_pos_init-repmat(mean(true_pos_init),3,1); %Transform the true initial positions, so the center is at (0,0)
rec_pos_trans_init=rec_pos_trans_init-repmat(mean(rec_pos_trans_init),3,1); %Center the reconstructed initial points

%Find rotation (and reflection) that matches the reconstructed and true locations of the reference points

%First do rough grid search of possible rotations and reflections
dthetas=linspace(0,7/4*pi,7);
sse_min=Inf;
fmin=Inf;
%Test all possible reflection combinations
for i=1:2
    temp_rec_pos_trans_init=rec_pos_trans_init;
    if i==2
        temp_rec_pos_trans_init(:,1)=-rec_pos_trans_init(:,1);
    end
    %Test all rotations given the current reflection
    for j=1:7
        dtheta=dthetas(j);
        sse=rotate_func_2d(dtheta,temp_rec_pos_trans_init,true_pos_trans_init); %Get the sum-squared error or the rotation
        if sse<sse_min %If this is the best rotation for this reflection, use these parameters as initial parameters for a more thorough search
            sse_min=sse;
            inits=dtheta;
        end     
    end
    %Find the best rotation (using fmincon instead of a gridsearch) for the
    %given reflection (using the gridsearch optimal parameters as the initial
    %parameters)
    [q,fval]=fmincon(@(x) rotate_func_2d(x,temp_rec_pos_trans_init,true_pos_trans_init),inits,[],[],[],[],0,2*pi);
    if fval<fmin
        fmin=fval; %best sum squared error
        minVal=i; %best reflection param
        minQ=q; %best rotation param
    end
end

%Based on the optimal values found above, do the reflection and rotation

%Do the optimal reflection
rec_pos_trans=rec_pos;
if minVal==2
    rec_pos_trans(:,1)=-rec_pos_trans(:,1);
end

%Do the optimal rotation
theta=minQ;
R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
rec_pos_trans=rec_pos_trans*R;


%% Test Reconstruction Accuracy (Using R-value)

%Compare pairwise distances between all points originally and in reconstruction
true_pos=[init_x init_y];
true_dists=pdist(true_pos);
model_dists=pdist(rec_pos_trans);
r2=corr(true_dists',model_dists');

%Note: You can easily also determine a mean error like in the other applications

%% Plot Initial and Reconstructed Points

if plottoggle
%Create colormap 
x=size(dists.images.original,1);
y=size(dists.images.original,2);
xVal=init_x/x;
yVal=init_y/y;
xVal=xVal(:);
yVal=yVal(:);
colors_xy = [.6*ones(size(xVal)),xVal,yVal];

%Plot Initial Points
init_x_trans=init_x-mean(init_x);
init_y_trans=init_y-mean(init_y);

figure; scatter(init_x_trans,init_y_trans,500,colors_xy,'.')
xlim([-700 700]); ylim([-700 700]);
%     set(gca,'xtick',[],'ytick',[]); box on;
set(gca,'visible','off')

%Plot Reconstructed Points
figure; scatter(rec_pos_trans(:,1),rec_pos_trans(:,2),500,colors_xy,'.')
% xlim([-700 700]); ylim([-700 700]);
%     set(gca,'xtick',[],'ytick',[]); box on;
set(gca,'visible','off')
end


%% Plot Initial Chemical Map

plate = dists.images.original;
if plottoggle == 1
    figure;
    imagesc(plate); colormap('gray');
    axis off;
end

%% Reconstruct Chemical Map

%First find the chemical compositions at the starting location of each
%pioneer cell (the chemical concentration it would have encoded)
%At the reconstructed locations of the pioneer cells, we thus know the
%chemical concentation.
%Then, we can fill in all the pixels of the chemical map using interpolation

%Create a vector of the chemical compositions of each cell 
%(taken from the image value at the starting pixels)
init_plate=zeros(1,length(rec_pos));
for i=1:length(rec_pos)
    prob=plate(init_y(i),init_x(i)); %This is the chemical concentration at the starting location of pioneer cell i.
    init_plate(i)=mean(rand(1,bclength)<prob); 
    %This (above line) is the chemical concentration that the pioneer cell would encode,
    %given the number of basepairs. Each basepair has a probability of
    %having a different nucleotide that is proportional to the chemical
    %concentration.                                                
end

%Create an empty chemical map where the pixels are in the range of the
%reconstructed locations
recreate_x=rec_pos_trans(:,1)-min(rec_pos_trans(:,1));
recreate_y=rec_pos_trans(:,2)-min(rec_pos_trans(:,2));
[chem_map_x,chem_map_y]=meshgrid(linspace(1,max(recreate_x),1000),linspace(1,max(recreate_y),1000));

%Interpolate chemical concentrations over all the pixels
interpmethod='linear';
fprintf('Interpolating data using %s method...\n',interpmethod); tic;
%At chemical map positions [chem_map_x, chem_map_y], interpolate the
%chemical concentration using the known concentrations, "init_plate", at
%the known locations [recreate_x,recreate_y]
chem_map=griddata(recreate_x,recreate_y,init_plate,chem_map_x(:),chem_map_y(:),interpmethod);
chem_map=reshape(chem_map,size(chem_map_x));

%% Plot Reconstructed Chemical Map

if plottoggle == 1
    figure;
    % subplot(1,2,2);
    imagesc(chem_map); colormap('gray'); %title(['Reconstructed Image, ' upper(interpmethod(1)) lower(interpmethod(2:end)) ' Interpolation']);
    axis off;
end
