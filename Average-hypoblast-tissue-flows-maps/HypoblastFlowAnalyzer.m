%% HypoblastFlowsAnalyzer
% Collect all hypoblast flows data (and their decomposition)
% Generate archetypal flow maps.

close all
clearvars

%% Parameters //////////////////////////////////////////////////////////////////////////////////////////////////////////
% Path where to find the data (one folder per animal, called ['Embryo' i])
Path='\\gaia.pasteur.fr\MVS_DATA1\Aurelien\PAPER\23-02-24_Hypoblast-paper\CodeAvailability\PIV-analysis_Archetypal-maps\Example';
% Path where to store the pooled data
PathOut='\\gaia.pasteur.fr\MVS_DATA1\Aurelien\PAPER\23-02-24_Hypoblast-paper\CodeAvailability\PIV-analysis_Archetypal-maps\Data';
% IDs of the embryo to analyze (here, only one movie is used as an example,
% but the real dataset contains 9 movies)
EmbryosIDs=1:1;

% Boundaries and intervals for the 3D interpolation
% Time (h)
tmin=2;
tmax=12;
tstep=0.5;
% X (in % of distance between the landmarks)
Xmin=-150;
Xmax=150;
Xstep=10;
% Y (in % of distance between the landmarks)
Ymin=-150;
Ymax=150;
Ystep=10;

%% Code ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
% Initialization of the pooling matrices
PoolX=nan(length(Xmin:Xstep:Xmax),length(Ymin:Ystep:Ymax),length(tmin:tstep:tmax),length(EmbryosIDs));
PoolY=nan(length(Xmin:Xstep:Xmax),length(Ymin:Ystep:Ymax),length(tmin:tstep:tmax),length(EmbryosIDs));
PoolDivX=nan(length(Xmin:Xstep:Xmax),length(Ymin:Ystep:Ymax),length(tmin:tstep:tmax),length(EmbryosIDs));
PoolDivY=nan(length(Xmin:Xstep:Xmax),length(Ymin:Ystep:Ymax),length(tmin:tstep:tmax),length(EmbryosIDs));
PoolRotatX=nan(length(Xmin:Xstep:Xmax),length(Ymin:Ystep:Ymax),length(tmin:tstep:tmax),length(EmbryosIDs));
PoolRotatY=nan(length(Xmin:Xstep:Xmax),length(Ymin:Ystep:Ymax),length(tmin:tstep:tmax),length(EmbryosIDs));
PoolLandmark100=nan(length(EmbryosIDs));

% For each embryo
for embryo=EmbryosIDs
    % Load data //////////////////////////////////////////////////////////
    % Load timing vector
    Timing=readmatrix([Path filesep 'Embryo' num2str(embryo) filesep 'Timing.csv']);
    % Load spatial landmarks coordinates
    SpatialLandmarks=readmatrix([Path filesep 'Embryo' num2str(embryo) filesep 'CoordinatesCentersRotation.csv']);
    % Load pixel size (in µm)
    Pixelsize=readmatrix([Path filesep 'Embryo' num2str(embryo) filesep 'Pixelsize.csv']);

    % Open the mask
    FileTif=[Path filesep 'Embryo' num2str(embryo) filesep 'Mask.tif'];
    InfoImage=imfinfo(FileTif);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    Mask=zeros(nImage,mImage,NumberImages,'uint16');
    TifLink = Tiff(FileTif, 'r');
    for j=1:NumberImages
        TifLink.setDirectory(j);
        Mask(:,:,j)=TifLink.read();
    end
    TifLink.close();

    % Deduce number of timepoints, height and width of the image
    Height=size(Mask,1);
    Width=size(Mask,2);
    NumberTimepoints=size(Mask,3);

    % Read PIV data
    Spacing=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(1) '/spacing']);
    xmin=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(1) '/xmin']);
    xmax=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(1) '/xmax']);
    ymin=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(1) '/ymin']);
    ymax=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(1) '/ymax']);

    % Initializations
    BinX=Spacing:Spacing:ymax-Spacing/2;
    BinY=Spacing:Spacing:xmax-Spacing/2;
    SpeedX=nan(length(BinY),length(BinX),NumberTimepoints-1);
    SpeedY=nan(length(BinY),length(BinX),NumberTimepoints-1);
    SpeedDivX=nan(length(BinY),length(BinX),NumberTimepoints-1);
    SpeedDivY=nan(length(BinY),length(BinX),NumberTimepoints-1);

    % Fill up the speed matrices frame after frame
    for t=0:NumberTimepoints-2
        % Load the right timing in .h5 file
        dx=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(t) '/dx']);
        dy=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(t) '/dy']);
        dx1=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(t) '/dx1']);
        dy1=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(t) '/dy1']);
        x=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(t) '/x']);
        y=h5read([Path filesep 'Embryo' num2str(embryo) filesep 'data' filesep 'MAX-flows.h5'],['/' num2str(t) '/y']);

        % Put the values in the appropiate place in the matrix
        for arrow=1:length(x)
            indexX=find(BinX==x(arrow));
            indexY=find(BinY==y(arrow));
            if ~isempty(indexX) && ~isempty(indexY)
                SpeedX(indexY,indexX,t+1)=dx(arrow);
                SpeedY(indexY,indexX,t+1)=dy(arrow);
                SpeedDivX(indexY,indexX,t+1)=dx1(arrow);
                SpeedDivY(indexY,indexX,t+1)=dy1(arrow);
            end
        end
    end

    % Convert speeds in µm/h
    SpeedX=SpeedX*Pixelsize/(Timing(2)-Timing(1));
    SpeedY=SpeedY*Pixelsize/(Timing(2)-Timing(1));
    SpeedDivX=SpeedDivX*Pixelsize/(Timing(2)-Timing(1));
    SpeedDivY=SpeedDivY*Pixelsize/(Timing(2)-Timing(1));


    % Apply mask /////////////////////////////////////////////////////////
    % Convert the mask so as to match time intervals
    Mask2=mean(cat(4,Mask(:,:,1:end-1),Mask(:,:,2:end)),4);
    Mask2(Mask2~=0)=1;
    % Sample the mask according to PIV boxes
    Mask3=zeros(length(BinY),length(BinX),size(Mask2,3));
    for t=1:size(Mask2,3)
        for y=1:length(BinY)
            for x=1:length(BinX)
                temp=Mask2(BinX(x)-Spacing/2:BinX(x)+Spacing/2,BinY(y)-Spacing/2:BinY(y)+Spacing/2,t);
                Covering=100*sum(temp(:)==1)/length(temp(:));
                if Covering>=50
                    Mask3(x,y,t)=1;
                end
            end
        end
    end
    % Apply mask to PIV data
    SpeedX(Mask3==0)=nan;
    SpeedY(Mask3==0)=nan;
    SpeedDivX(Mask3==0)=nan;
    SpeedDivY(Mask3==0)=nan;


    % Spatio-temporal alignment //////////////////////////////////////////
    % Convert the landmarks to real image size
    SpatialLandmarks=SpatialLandmarks/600*max([Height,Width]);
    % Deduce timing of the PIV intervals from frames timing
    TimingPIV=Timing+(Timing(2)-Timing(1));
    TimingPIV(end)=[];
    % Convert coordinates of the PIV box in coordinates based on landmarks
    % (left landmerk=[-50,0], right landmark=[50,0])
    Landmark100=SpatialLandmarks(3)-SpatialLandmarks(1);
    LandmarkY=mean([SpatialLandmarks(2),SpatialLandmarks(4)]);
    BinY=100*(BinY-LandmarkY)/Landmark100;
    BinX=100*(BinX-SpatialLandmarks(1)-Landmark100/2)/Landmark100;
    % 3D interpolation
    [X,Y,Z]=meshgrid(BinX,BinY,TimingPIV);
    [Xq,Yq,Zq]=meshgrid(Xmin:Xstep:Xmax,Ymin:Ystep:Ymax,tmin:tstep:tmax);
    InterpSpeedX=interp3(X,Y,Z,SpeedX,Xq,Yq,Zq);
    InterpSpeedY=interp3(X,Y,Z,SpeedY,Xq,Yq,Zq);
    InterpSpeedDivX=interp3(X,Y,Z,SpeedDivX,Xq,Yq,Zq);
    InterpSpeedDivY=interp3(X,Y,Z,SpeedDivY,Xq,Yq,Zq);

    % Pooling all animals and normalization by the distance between
    % landmarks
    Landmark100=Landmark100*Pixelsize;
    PoolX(:,:,:,embryo)=InterpSpeedX*100/Landmark100;
    PoolY(:,:,:,embryo)=InterpSpeedY*100/Landmark100;
    PoolDivX(:,:,:,embryo)=InterpSpeedDivX*100/Landmark100;
    PoolDivY(:,:,:,embryo)=InterpSpeedDivY*100/Landmark100;
    PoolRotatX(:,:,:,embryo)=InterpSpeedX*100/Landmark100-InterpSpeedDivX*100/Landmark100;
    PoolRotatY(:,:,:,embryo)=InterpSpeedY*100/Landmark100-InterpSpeedDivY*100/Landmark100;
    PoolLandmark100(embryo)=Landmark100;

    disp(['Embryo' num2str(embryo) ' successfully processed!']);
end

% % Save pooled data (here, this section is commented as only one example
% is provided, but the output for the 9 movies is provided in the "Data"
% folder)
% save([PathOut filesep 'PoolX.mat'],'PoolX');
% save([PathOut filesep 'PoolY.mat'],'PoolY');
% save([PathOut filesep 'PoolDivX.mat'],'PoolDivX');
% save([PathOut filesep 'PoolDivY.mat'],'PoolDivY');
% save([PathOut filesep 'PoolRotatX.mat'],'PoolRotatX');
% save([PathOut filesep 'PoolRotatY.mat'],'PoolRotatY');


%% Generate average maps

% Load data
load([PathOut filesep 'PoolX.mat']);
load([PathOut filesep 'PoolY.mat']);
load([PathOut filesep 'PoolDivX.mat']);
load([PathOut filesep 'PoolDivY.mat']);
load([PathOut filesep 'PoolRotatX.mat']);
load([PathOut filesep 'PoolRotatY.mat']);

% Generate a covering filter
Covering=sum(~isnan(PoolX),4);
VectT=tmin:tstep:tmax;
T4h=find(VectT==4);
T8h=find(VectT==8);
T12h=find(VectT==12);

% Average of all animals for each time point
AverageX=mean(PoolX,4,'omitnan');
AverageY=mean(PoolY,4,'omitnan');
AverageX(Covering<=3)=nan;
AverageY(Covering<=3)=nan;
AverageRotatX=mean(PoolRotatX,4,'omitnan');
AverageRotatY=mean(PoolRotatY,4,'omitnan');
AverageRotatX(Covering<=3)=nan;
AverageRotatY(Covering<=3)=nan;
AverageDivX=mean(PoolDivX,4,'omitnan');
AverageDivY=mean(PoolDivY,4,'omitnan');
AverageDivX(Covering<=3)=nan;
AverageDivY(Covering<=3)=nan;

%% Calculate vorticity and divergence at 8h and 12h

% Create a color bar orange to magenta
mycolormap = customcolormap([0 0.5 1], [255/255 191/255 0; 1 1 1; 1 0 1]);

% 8h
step=2;
T=8;
t=find(VectT==T);
Time=VectT(t);
VX=nanmean(AverageX(:,:,t-step:t+step),3);
VY=nanmean(AverageY(:,:,t-step:t+step),3);
% Vorticity and divergence calculation
Vorticity=nan(size(VX));
Divergence=nan(size(VX));
for x=2:size(VX,1)-1
    for y=2:size(VX,2)-1
        if ~isnan(VX(x,y)) & ~isnan(VX(x-1,y)) & ~isnan(VX(x,y-1)) & ~isnan(VX(x+1,y)) & ~isnan(VX(x,y+1))
            % Calculation of vorticity and divergence
            Vorticity(x,y)=( VY(x,y-1)-VY(x,y+1) ) / (2*Xstep) - ( VX(x-1,y) - VX(x+1,y) ) / (2*Xstep);
            Divergence(x,y)=( VX(x,y-1)-VX(x,y+1) ) / (2*Xstep) + ( VY(x-1,y)-VY(x+1,y) ) / (2*Xstep);
        end
    end
end
hfig=figure(1);
Vorticity(isnan(Vorticity))=0;
imagesc(Vorticity)
caxis([-1 1])
colormap(mycolormap)
daspect([1 1 1])
xticklabels([]);
yticklabels([]);
colorbar
saveas(hfig,[Path filesep '7-9h(vorticity).pdf']);
hfig=figure(2);
Divergence(isnan(Divergence))=0;
imagesc(Divergence)
caxis([-0.5 0.5])
colormap(redblue)
daspect([1 1 1])
xticklabels([]);
yticklabels([]);
colorbar
saveas(hfig,[PathOut filesep '7-9h(divergence).pdf']);

% 11h
step=2;
T=11;
t=find(VectT==T);
Time=VectT(t);
VX=nanmean(AverageX(:,:,t-step:t+step),3);
VY=nanmean(AverageY(:,:,t-step:t+step),3);
% Vorticity and divergence calculation
Vorticity=nan(size(VX));
Divergence=nan(size(VX));
for x=2:size(VX,1)-1
    for y=2:size(VX,2)-1
        if ~isnan(VX(x,y)) & ~isnan(VX(x-1,y)) & ~isnan(VX(x,y-1)) & ~isnan(VX(x+1,y)) & ~isnan(VX(x,y+1))
            % Calculation of vorticity and divergence
            Vorticity(x,y)=( VY(x,y-1)-VY(x,y+1) ) / (2*Xstep) - ( VX(x-1,y) - VX(x+1,y) ) / (2*Xstep);
            Divergence(x,y)=( VX(x,y-1)-VX(x,y+1) ) / (2*Xstep) + ( VY(x-1,y)-VY(x+1,y) ) / (2*Xstep);
        end
    end
end
hfig=figure(3);
Vorticity(isnan(Vorticity))=0;
imagesc(Vorticity)
caxis([-1 1])
colormap(mycolormap)
daspect([1 1 1])
xticklabels([]);
yticklabels([]);
saveas(hfig,[Path filesep '10-12h(vorticity).pdf']);
hfig=figure(4);
Divergence(isnan(Divergence))=0;
imagesc(Divergence)
caxis([-0.5 0.5])
colormap(redblue)
daspect([1 1 1])
xticklabels([]);
yticklabels([]);
saveas(hfig,[PathOut filesep '10-12h(divergence).pdf']);


% Save 7-9h and 10-12h vector fields for overlaying on Vorticity and
% Divergence
t=find(VectT==8);
Time=VectT(t);
factor=2;
hfig=figure;
subplot(1,3,1)
[X,Y]=meshgrid(Xmin:Xstep:Xmax,Ymin:Ystep:Ymax);
quiver(X,Y,factor*nanmean(AverageX(:,:,t-step:t+step),3),factor*nanmean(AverageY(:,:,t-step:t+step),3),'AutoScale','off','ShowArrowHead','off');
set(gca, 'YDir', 'reverse');
daspect([1 1 1])
title([num2str(Time) 'h'])
xlim([-150 150])
ylim([-150 150])
xticklabels([]);
yticklabels([]);
subplot(1,3,2)
quiver(X,Y,factor*nanmean(AverageRotatX(:,:,t-step:t+step),3),factor*nanmean(AverageRotatY(:,:,t-step:t+step),3),'AutoScale','off','ShowArrowHead','off');
set(gca, 'YDir', 'reverse');
daspect([1 1 1])
title([num2str(Time) 'h'])
xlim([-150 150])
ylim([-150 150])
xticklabels([]);
yticklabels([]);
subplot(1,3,3)
quiver(X,Y,factor*nanmean(AverageDivX(:,:,t-step:t+step),3),factor*nanmean(AverageDivY(:,:,t-step:t+step),3),'AutoScale','off','ShowArrowHead','off');
set(gca, 'YDir', 'reverse');
daspect([1 1 1])
title([num2str(Time) 'h'])
xlim([-150 150])
ylim([-150 150])
xticklabels([]);
yticklabels([]);
set(gcf,'Position',[50, 50, 600,200]);
saveas(hfig,[PathOut filesep '7-9h.pdf']);

t=find(VectT==11);
Time=VectT(t);
hfig=figure;
subplot(1,3,1)
[X,Y]=meshgrid(Xmin:Xstep:Xmax,Ymin:Ystep:Ymax);
quiver(X,Y,factor*nanmean(AverageX(:,:,t-step:t+step),3),factor*nanmean(AverageY(:,:,t-step:t+step),3),'AutoScale','off','ShowArrowHead','off');
set(gca, 'YDir', 'reverse');
daspect([1 1 1])
title([num2str(Time) 'h'])
xlim([-150 150])
ylim([-150 150])
xticklabels([]);
yticklabels([]);
subplot(1,3,2)
quiver(X,Y,factor*nanmean(AverageRotatX(:,:,t-step:t+step),3),factor*nanmean(AverageRotatY(:,:,t-step:t+step),3),'AutoScale','off','ShowArrowHead','off');
set(gca, 'YDir', 'reverse');
daspect([1 1 1])
title([num2str(Time) 'h'])
xlim([-150 150])
ylim([-150 150])
xticklabels([]);
yticklabels([]);
subplot(1,3,3)
quiver(X,Y,factor*nanmean(AverageDivX(:,:,t-step:t+step),3),factor*nanmean(AverageDivY(:,:,t-step:t+step),3),'AutoScale','off','ShowArrowHead','off');
set(gca, 'YDir', 'reverse');
daspect([1 1 1])
title([num2str(Time) 'h'])
xlim([-150 150])
ylim([-150 150])
xticklabels([]);
yticklabels([]);
set(gcf,'Position',[50, 50, 600,200]);
saveas(hfig,[PathOut filesep '10-12h.pdf']);