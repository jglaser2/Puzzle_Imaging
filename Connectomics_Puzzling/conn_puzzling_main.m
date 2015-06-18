%%This script will do all of the connectomics puzzling applications from
%%the manuscript. Its sections are:
%Layer 3: Data for Fig 4C and D, and Plot for Fig 4A
%Layers 2 and 3: Data for Fig 4C and D, and Plot for Fig 4B
%Changing sigma parameter: Data for Fig 4G
%Changing A parameter: Data for Fig 4F



%% Layer 3

%IF PLOT_FIG=0, THEN THIS CALCULATES THE RECONSTRUCTION ACCURACY FOR LAYER 3(DATA FOR FIGS 4C AND D)
%IF PLOT_FIG=1, THEN THIS PLOTS AN EXAMPLE RECONSTRUCTION OF LAYER 3 (FIG 4A)
plot_fig=0; 


%GET DATA FOR FIGS 4C,D FOR LAYER 3
if ~plot_fig

reps=3; %200

%Layer 3-->3 params (see manuscript)
x0_33=5.705;
p0_33=0.0720;
w_33=296.7;
A_33=270;
p33 = @(x) p0_33+A_33/w_33/sqrt(pi/2)*exp(-2*(x-x0_33).^2/w_33^2);

r2_vec=zeros(reps,1); %Vector keeping track of r-values between true and reconstructed pairwise distances between neurons
err_mean_vec=zeros(reps,1); %Vector keeping track of mean distance errors of the reconstruction
removed=zeros(reps,1); %Vector keeping track of how many neurons needed to be excluded because they weren't part of the main connected component

x=20; %Number of neurons that fit across simulation in x-dimension
y=20; %Number of neurons that fit across simulation in y-dimension
z=20; %Number of neurons that fit across simulation in z-dimension
scale=20; %distance between neurons

for rep=1:reps
    
    [ r2_vec(rep), err_mean_vec(rep), removed(rep) ] = conn_puzzling_func( x,y,z,scale,0,p33 );
    
end

% save('Fig4C_layer3_diff','r2_vec','err_vec','err_mean_vec','removed')
end


%PLOTS
if plot_fig
    conn_puzzling_func( x,y,z,scale,1,p33 );
end

%% Layers 2 and 3

%IF PLOT_FIG=0, THEN THIS CALCULATES THE RECONSTRUCTION ACCURACY FOR LAYERS 2 AND 3(DATA FOR FIGS 4C AND D)
%IF PLOT_FIG=1, THEN THIS PLOTS AN EXAMPLE RECONSTRUCTION OF LAYERS 2 AND 3 (FIG 4B)
plot_fig=0; 


%GET DATA FOR FIGS 4C,D FOR LAYERS 2 AND 3
if ~plot_fig
    
reps=3; %200


%Layer 3-->3 params
x0_33=5.705;
p0_33=0.0720;
w_33=296.7;
A_33=270;
p33 = @(x) p0_33+A_33/w_33/sqrt(pi/2)*exp(-2*(x-x0_33).^2/w_33^2);

%Layer2-->2 params
x0_22=-35.33;
p0_22=0.0714;
w_22=489.9;
A_22=459.3;
p22 = @(x) p0_22+A_22/w_22/sqrt(pi/2)*exp(-2*(x-x0_22).^2/w_22^2);

%Layer 2-->3 params
x0_23=-259.4;
p0_23=-0.0231;
w_23=620.8;
A_23=987.1;
p23 = @(x) p0_23+A_23/w_23/sqrt(pi/2)*exp(-2*(x-x0_23).^2/w_23^2);

r2_vec=zeros(reps,1); %Vector keeping track of r-values between true and reconstructed pairwise distances between neurons
err_mean_vec=zeros(reps,1); %Vector keeping track of mean distance errors of the reconstruction
removed=zeros(reps,1); %Vector keeping track of how many neurons needed to be excluded because they weren't part of the main connected component

x=20; %Number of neurons that fit across simulation in x-dimension
y=20; %Number of neurons that fit across simulation in y-dimension
z=20; %Number of neurons that fit across simulation in z-dimension
scale=20; %distance between neurons

for rep=1:reps
    
    [ r2_vec(rep), err_mean_vec(rep), removed(rep) ] = conn_puzzling_func( x,y,z,0,scale,p33,p22,p23 );
    
end

% save('Fig4C_layer23_diff','r2_vec','err_vec','err_mean_vec','removed')

end


%PLOTS
if plot_fig
    conn_puzzling_func( x,y,z,scale,1,p33,p22,p23 );
end

%% Changing sigma parameter

%Number of simulations for each sigma value
reps=10;

%Changing sigma parameter
sigmas=10*(5:5:200);

%Initializations
r2=zeros(length(sigmas),reps); %Matrix keeping track of r-values between true and reconstructed pairwise distances between neurons for each sigma
err_mean=zeros(length(sigmas),reps); % Matrix keeping track of mean distance errors of the reconstruction for each sigma
removed=r2; %Matrix keeping track of how many neurons needed to be excluded because they weren't part of the main connected component, for each sigma

x=20; %Number of neurons that fit across simulation in x-dimension
y=20; %Number of neurons that fit across simulation in y-dimension
z=20; %Number of neurons that fit across simulation in z-dimension
scale=20; %distance between neurons


tic
for jj=1:length(sigmas)
    
    %Probability function parameters
    w=2*sigmas(jj); %width is 2 sigma
    x0=5.705;
    p0=0.0720;
    A=270/296.7*w;    %A/w is constant
    p = @(x) p0+A/w/sqrt(pi/2)*exp(-2*(x-x0).^2/w^2);
    
    
    for kk=1:reps
        
        [ r2(jj,kk), err_mean(jj,kk), removed(jj,kk) ] = conn_puzzling_func( x,y,z,scale,0,p);
        
    end
end

% save('Fig4_sigma_diff','r2','err_mean','removed')

%% Changing baseline probability (A) parameter

%Number of simulations for each sigma value
reps=10;

%Changing baseline parameter
base=.05:.05:1;

%Other probability function parameters
x0=5.705;
p0=0.0720;
w=296.7;
As=w*sqrt(pi/2)*(base-p0);

%Initializations
r2=zeros(length(base),reps); %Matrix keeping track of r-values between true and reconstructed pairwise distances between neurons for each sigma
err_mean=zeros(length(base),reps); % Matrix keeping track of mean distance errors of the reconstruction for each sigma
removed=r2; %Matrix keeping track of how many neurons needed to be excluded because they weren't part of the main connected component, for each sigma

x=20; %Number of neurons that fit across simulation in x-dimension
y=20; %Number of neurons that fit across simulation in y-dimension
z=20; %Number of neurons that fit across simulation in z-dimension
scale=20; %distance between neurons


tic
for jj=1:length(base)
    
    %Probability function
    A=As(jj);
    p = @(x) p0+A/w/sqrt(pi/2)*exp(-2*(x-x0).^2/w^2);
    
    for kk=1:reps
        
        [ r2(jj,kk), err_mean(jj,kk), removed(jj,kk) ] = conn_puzzling_func( x,y,z,scale,0,p);
        
    end
end

% save('Fig4_Aparam_diff','r2','err_mean','rl')