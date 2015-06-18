function dists = chem_sim(pion_cells,probconj)
%This function does the chemical puzzling simulation that the results in the paper come from.
%The paper also mentions a more complex (but more timely) simulation.
%Please email to request this code (although it's currently not commented well...)


%As input, it takes "pion_cells", which can either be the number of pioneer
%cells or the density of pioneer cells. It also takes "probconj", which is
%the probability of conjugation between an F+ and F- cell.

%As output, it gives a struct, "dists" which gives the locations of the
%pioneer cells (dists.init), the image/plate being used (dists.images), and
%the connectivity matrix between pioneer cells (dists.josh)


%Make default conjugation probability 100%
if nargin<2
    probconj=1;
end


%% Load Image/Plate

%Load the grayscale image (that represents the chemical concentrations across the plate)
imdata = double(imread('P3.jpg'));
imdata=imdata/max(imdata(:)); %Scale to have maximum of 1

%Get information about the size of the plate
num_x = size(imdata,1);
num_y = size(imdata,2);
numpix=num_x*num_y;

%% Other basics

%If pion_cells is entered as a cell density, convert it to the number of pioneer cells
if pion_cells < 1
    pion_cells = uint64(floor(numpix*pion_cells));
end

%Whether each pioneer cell is F+ or F- (conj=1 or 0)
conj = randi(2,pion_cells,1)-1;

%% Generate the unique initial locations for all the pioneer cells

% These are the known (reference) pioneer cell positions (x,y) as fraction of entire plate
inits = [0.25 0.25; 0.25 0.75; 0.75 0.25];

%Generate a list of many unique possible x and y positions 
%(there should be more than enough unique combinations, but there's a while loop just to make sure).
uM = [];
while(size(uM,1)<pion_cells-size(inits,1))
    initx = randi(num_x,10*pion_cells,1); 
    inity = randi(num_y,10*pion_cells,1);
    uM = unique([initx inity],'rows');
end

%Select the unique location for all the pioneer cells (from the above too-long list)
uM = sortrows(uM(randperm(size(uM,1),pion_cells-size(inits,1)),:));
% reference cells are placed at predefined locations (see top of script)
init_pioneers(:,1) = floor(num_x*inits(:,1));
init_pioneers(:,2) = floor(num_y*inits(:,2));
uM = [init_pioneers; uM];
initx = uM(:,1);
inity = uM(:,2);

%% Find the "connections" between pioneer cells (conjugations between colonies)

%First, assuming each pioneer cell's colony grows outward at an equal rate, find
%what pioneer cell's colony will end up occuping each "pixel" of the plate.
%The colony that gets to a pixel first is assumed to occupy that pixel.
%We do this by making a 2d Gaussian surrounding each pioneer cell's initial
%location. Then for each pixel, we find the Gaussian that has the largest
%value of its pdf.

[x0,y0]=meshgrid(1:num_x, 1:num_x);
pos=[x0(:) y0(:)];

MU=[initx(1) inity(1)];
SIGMA=[10,10];
y = mvnpdf(pos,MU,SIGMA);
maxy=reshape(y,size(y0)); %This keeps track of the maximum value of any Gaussian at each pixel
GaussMaxMat=ones(num_x,num_y); %This is the matrix that keeps track of which pioneer cell's colony ends up in each pixel

for i=2:pion_cells %Loop through pioneer cells
    MU=[initx(i) inity(i)];    
    y = mvnpdf(pos,MU,SIGMA);
    y=reshape(y,size(y0));
    [a]=find(y>maxy); %Find the pixels that this pioneer cell's Gaussian has a larger value in compared to previously tested pioneer cells
    maxy(a)=y(a); %Update the maximums
    GaussMaxMat(a)=i; %Update which pioneer cell ends up in each pixel (given the pioneer cells we have looped through)
end

%Next, find which pioneer cell colonies conjugate with each other

%To do this, we first find all "border" pixels, meaning that there are different
%pioneer cell colonies in adjacent pixels.

diffMatV=(GaussMaxMat((1:end-1),:)-GaussMaxMat((2:end),:))~=0; %Pixels that have different colonies vertically (above or below)
diffMatH=(GaussMaxMat(:,(1:end-1))-GaussMaxMat(:,(2:end)))~=0; %Pixels that have different colonies horizontally (left or right)

%Then, for each border pixel, give a connection between the occupiying
%pioneer cell colonies with probability "probconj".

conn=zeros(pion_cells); %Matrix of connections

%Determine connections for vertical border pixels
[a,b,~]=find(diffMatV);
for i=1:length(a)
    c1=GaussMaxMat(a(i),b(i)); %Find which pioneer cell colonies are bordering
    c2=GaussMaxMat(a(i)+1,b(i));
    if rand<probconj
        conn(c1,c2)=conn(c1,c2)+1;
        conn(c2,c1)=conn(c2,c1)+1;
    end
end

%Determine connections for horizontal border pixels
[a,b,~]=find(diffMatH);
for i=1:length(a)
    c1=GaussMaxMat(a(i),b(i)+1);
    c2=GaussMaxMat(a(i),b(i));
    if rand<probconj
        conn(c1,c2)=conn(c1,c2)+1;
        conn(c2,c1)=conn(c2,c1)+1;
    end
end

%Finally, there can only be connections between colonies if one is an F+
%and another is an F-. So, remove all connections between F+ and F+ or
%between F- and F-

conj_xor=squareform(pdist(conj)); %Determines whether every two pioneer cells are conjugation compatible (one is F+ and the other is F-).
conn2=conn.*conj_xor; %Remove the incompatible conjugations

%% Save output

dists.images.original = imdata;
dists.init.x = initx;
dists.init.y = inity;
% dists.conj_xor=conj_xor;
% dists.conn=conn;
dists.josh = conn2;