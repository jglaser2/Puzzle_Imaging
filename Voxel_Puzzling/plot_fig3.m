%Script to make results plots in Figure 3
%Note that "diff" refers to SDM and "Tri" refers to ULI (names changed since the beginning of the project)


%% Fig 3 C,D,E
voxsize=2:15;
voxsize2=2:10;

load('Fig2_diff_param.mat')
r2_lap=r2;
err_lap=err_mean;
rem_lap=removed;
load('Fig2_tri_param_5.mat')
r2_tri=r2;
err_tri=err_mean;
rem_tri=removed;

mean_r2_tri=zeros(size(r2_tri,1),1);
std_r2_tri=zeros(size(r2_tri,1),1);
mean_err_tri=zeros(size(r2_tri,1),1);
std_err_tri=zeros(size(r2_tri,1),1);
mean_r2_lap=zeros(size(r2_lap,1),1);
std_r2_lap=zeros(size(r2_lap,1),1);
mean_err_lap=zeros(size(r2_lap,1),1);
std_err_lap=zeros(size(r2_lap,1),1);
for i=1:size(r2_lap,1)
    idx=~isnan(r2_lap(i,:));
    mean_r2_lap(i)=mean(r2_lap(i,idx));
    std_r2_lap(i)=std(r2_lap(i,idx));
    mean_err_lap(i)=mean(err_lap(i,idx));
    std_err_lap(i)=std(err_lap(i,idx));
end
for i=1:size(r2_tri,1)
    idx=~isnan(r2_tri(i,:));
    mean_r2_tri(i)=mean(r2_tri(i,idx));
    std_r2_tri(i)=std(r2_tri(i,idx));
    mean_err_tri(i)=mean(err_tri(i,idx));
    std_err_tri(i)=std(err_tri(i,idx));
end

figure; errorbar(voxsize,mean_r2_lap,std_r2_lap)
hold on;
errorbar(voxsize2,mean_r2_tri,std_r2_tri,'r')
ylim([0 1])
xlim([0 10.2])
box off;

figure; errorbar(voxsize,mean_err_lap,std_err_lap)
hold on;
errorbar(voxsize2,mean_err_tri,std_err_tri,'r')
ylim([0 5])
xlim([0 10.2])
box off;

figure; errorbar(voxsize,voxsize.*mean_err_lap',voxsize.*std_err_lap')
hold on;
errorbar(voxsize2,voxsize2.*mean_err_tri',voxsize2.*std_err_tri','r')
xlim([0 10.2])
ylim([0 50])
box off;

% rem_total=[rem_lap rem_tri];
% mean_rem=mean(rem_total,2);
% perc_rem=mean_rem/8000;



%% Fig 3 F,G
% percremoved=0:.05:.75;
percremoved2=0:.05:.9;

load('Fig2_diff_holes.mat')
r2_lap=r2;
err_lap=err_mean;
rem_lap=removed;

load('Fig2_diff_holes_lots.mat')
r2_lap=[r2_lap; r2];
err_lap=[err_lap; err_mean];
rem_lap=[rem_lap; removed];


load('Fig2_tri_holes_lots.mat')
r2_tri=r2;
err_tri=err_mean;
rem_tri=removed;

for i=1:size(r2_tri,1)
    idx=~isnan(r2_tri(i,:));
    mean_r2_tri(i)=mean(r2_tri(i,idx));
    std_r2_tri(i)=std(r2_tri(i,idx));
    mean_err_tri(i)=mean(err_tri(i,idx));
    std_err_tri(i)=std(err_tri(i,idx));
end




figure; errorbar(percremoved2,mean(r2_lap'),std(r2_lap'))
hold on;
errorbar(percremoved2,mean_r2_tri,std_r2_tri,'r')
ylim([0 1])
xlim([-.01 0.91])
box off;

figure; errorbar(percremoved2,mean(err_lap'),std(err_lap'))
hold on;
errorbar(percremoved2,mean_err_tri,std_err_tri,'r')
xlim([-.01 0.91])
ylim([0 10])
box off;








