function [dist_vec]=dist_from_pt2(S,pt,numSteps)

%This function outputs dist_vec, the geodesic distance of all points to a given landmark point
%The inputs are "S" (the similarity matrix), "pt" (the index of the given
%landmark point, and "numSteps" (an estimate of the maximum geodesic
%distance to any point).

n=length(S); %number of total points

%reached is a vector that gives the points reached at the current step
%the initial condition (step 0) is that it has just reached the landmark point 
reached=zeros(1,n); 
reached(pt)=1;
reached=sparse(reached);

%reachedMat is a matrix that keeps track of which points have been reached
%after each number of steps.
reachedMat=zeros(numSteps,n); %numSteps is a guess to initialize the size of the matrix. It will grow in the following loop if the guess is too small.


for step=1:100000 %go through steps (100000 is an arbitrary large number - the loop will be broken out of before this)
    reached_old=reached; %reached_old gives the points reached on the previous step
    reached=reached_old*S; %Find what points have been reached at current step
    reachedMat(step,:)=reached'; %Fill in reachedMat for current step
    %If the points reached on the previous step are the same as those on
    %the current step, that means all the points have been reached. If so,
    %break out of the loop
    if(all((reached_old>0)==(reached>0))) 
        break
    end
end


%To find the geodesic distance of a point i to the landmark point, we need to
%find the first step at which reachedMat(:,i) is nonzero.
%To do this, we find the minimum nonzero values of each column of reachedMat 
%dist_vec gives the geodesic distance of all points to a given landmark point

reachedMat(~reachedMat)=Inf; %Make the zero values = infinity so we can get the minimum of nonzero entries in the next step
[~,dist_vec]=min(reachedMat);