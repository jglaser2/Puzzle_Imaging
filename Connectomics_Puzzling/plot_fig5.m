%This script generates the original plots used for Figure 5

%% Fig 5C 

load('Fig4C_layer3_diff.mat');
r2_vec3=r2_vec;
load('Fig4C_layer23_diff.mat');
r2_vec23=r2_vec;

r2mat=[r2_vec3 r2_vec23];

figure; boxplot(r2mat,'outliersize',0)
ylim([0.9 1])
box off;
% set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[])


%% Fig 5D

load('Fig4C_layer3_diff.mat');
err_mean_vec3=err_mean_vec;
load('Fig4C_layer23_diff.mat');
err_mean_vec23=err_mean_vec;

errmat=[err_mean_vec3 err_mean_vec23];

figure; boxplot(errmat,'outliersize',0)
ylim([0 50])
box off;
% set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[])

%% Fig 5E (Probability of connection plots)

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


figure; plot(p33(0:700))
hold on;
plot(p22(0:700),'r')
plot(p23(0:700),'g')
xlim([0 700])
ylim([0 1])
box off;
% set(gca,'XTickLabel',[],'YTickLabel',[])

%% Fig 5F

load('Fig4_Aparam_diff')


base=.05:.05:1;
mr2=mean(r2,2);
sr2=std(r2,0,2);
figure;
errorbar(base,mr2,sr2);
xlim([0 1.01])
ylim([0 1]);
% ylabel('R2 values')
% xlabel('Connection probability at distance 0')
% xlabel('Connection probability at distance 0')
% box off;
hold on;
mm=mean(err_mean,2);
sm=std(err_mean,0,2);

% figure;
errorbar(base,mm/300,sm/300,'r');
% xlim([0 1.01])
% ylim([0 300]);
% set(gca,'YAxisLocation','right')
% ylabel('Mean errors')
xlabel('Connection probability at distance 0')
% set(gca,'Visible','off')
% set(gca,'XTick',[],'YTick',[])
box off;


%% Fig 5G

load('Fig4_sigma_diff');

sigmas=10*(5:5:200);

mr2=mean(r2,2);
sr2=std(r2,0,2);
figure;
errorbar(sigmas,mr2,sr2);
xlim([0 2010])
ylim([0 1]);
% ylabel('R2 values')
% xlabel('Standard deviation of connection probability distribution (um)')
% box off;
hold on;

mm=mean(err_mean,2);
sm=std(err_mean,0,2);
% figure;
errorbar(sigmas,mm/300,sm/300,'r');
% xlim([0 2010])
% ylim([0 300]);
% set(gca,'YAxisLocation','right')
% ylabel('Mean errors')
% xlabel('Standard deviation of connection probability distribution (um)')
set(gca,'XTick',[],'YTick',[])
box off;
