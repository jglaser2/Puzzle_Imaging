function [rec_pos,pts]=uli_func(S,num_landmarks,dim)
%This function does the ULI (unweighted landmark isomap) algorithm to reduce dimensions
%The inputs are:
%S- the similarity matrix
%num_landmarks- the number of landmark points to use in the algorithm
%dim- how many dimensions the reconstructed positions are in (2 for
%chemical puzzling simulations, 3 for voxel puzzling)

%% Find geodesic distance of all things to landmark points

n=length(S); %number of points
numSteps=40; %Estimate of largest distance to a landmark point (if it's not accurate, it will just take a little longer)
Distmat=zeros(n,num_landmarks); % distance matrix from basis points
pts=zeros(1,num_landmarks); %keeps track of the indices of the landmark points

%First landmark point
pts(1)=randi(n); %Choose first landmark point randomly
b=dist_from_pt2(S,pts(1),numSteps); %Find the distance from all points to the landmark point
Distmat(:,1)=b; %Fill in the distance  matrix
Distmat(pts(1),1)=0; %Make the distance to itself be 0

%All other landmark points
for i=2:num_landmarks
    %Select landmark point using the maxmin algorithm (see manuscript)
    a=Distmat(:,1:i-1);
    mins=min(a,[],2);
    [~,idx]=max(mins);
    pts(i)=idx;
    b=dist_from_pt2(S,pts(i),numSteps); %Find the distance from all points to the landmark point
    Distmat(:,i)=b; %Fill in the distance  matrix
    Distmat(pts(i),i)=0; %Make the distance to itself be 0
end

%% Place the landmark points using MDS

rec_pos=zeros(n,dim); %Matrix of the recovered positions

D=Distmat(pts,:); %Geodesic distance from all points to the landmark points (written as gamma in the manuscript)

%Do MDS on the landmark points 
L=cmdscale(D,dim); %L from the manuscript

%% Triangulate all points based on the distances to the landmark points

D2=D.^2; %Delta^l from the manuscript
Distmat2=Distmat.^2; %Delta from the manuscript
distmean=mean(D2); %mean(Delta^l) from the manuscript
L_sharp=pinv(L); %L^sharp from the manuscript

for i=1:n
    if any(i==pts) %If the point is a landmark point, it's recovered position is that determined by MDS in the previous step
        rec_pos(i,:)=L(find(i==pts),:);        
    else
        distvec=Distmat2(i,:); %Delta_a in the manuscript (although here it's row i rather than row a)
        rec_pos(i,:)=.5*L_sharp*transpose((distmean-distvec)); %Triangulate!
    end   
end

%% 

%Below is code for doing MDS by the method written out in the paper instead
%of using cmdscale (and determining L_sharp is the way described in the
%paper)

% distmean=mean(D2);
% distmeanmean=mean(distmean);
% B=zeros(3,3);
% for i=1:3
%     for j=1:3
%         B(i,j)=-(D2(i,j)-distmean(i)-distmean(j)+distmeanmean)/2;
%     end
% end
% [v,e]=eigs(B,2);
% L2=v*sqrt(e);
%
% L_sharp=transpose(v/sqrt(e));
% distmean=mean(D2);
% distvec=Distmat2(i,:);
% rec_pos=.5*L_sharp*transpose((distmean-distvec));
