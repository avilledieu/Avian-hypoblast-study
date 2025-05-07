%% CompareHypoblastEpiblastFlows
% Upload manual hypoblast tracking and compare them to epiblast
% PIV-calculated flows
clearvars
close all

%% Parameters /////////////////////////////////////////////////////////////
% Path where to find subfolders containing chimera data ("Example" folder)
Path='\\gaia.pasteur.fr\MVS_DATA1\Aurelien\PAPER\23-02-24_Hypoblast-paper\CodeAvailability\PIV-analysis_Compare-with-hypoblast-flows\Example';
% Path where to find pooled data ("Data" folder)
PathData='\\gaia.pasteur.fr\MVS_DATA1\Aurelien\PAPER\23-02-24_Hypoblast-paper\CodeAvailability\PIV-analysis_Compare-with-hypoblast-flows\Data';
% Indexes of chimeras to treat (here, 1 example is provided, but the total
% dataset contains 10 animals)
Indexes=1:1;

% Smoothing radius for manual tracking smoothing (in h)
SmoothingRadiusT=0.5;
% /////////////////////////////////////////////////////////////////////////

%% Code ///////////////////////////////////////////////////////////////////
% Initialization of vectors storing the average displacement between
% hypoblast and epiblast and average displacement of hypoblast for each
% animal
PoolAverageDiff=nan(4,length(Indexes));
% Initialization of vectors pooling aligned hypoblast and difference of
% displacements
PoolX=[];
PoolY=[];
PoolAnimal=[];
PoolHypoX=[];
PoolHypoY=[];
PoolDiffX=[];
PoolDiffY=[];
PoolEpiX=[];
PoolEpiY=[];

% Treat each chimera after another/////////////////////////////////////////
for animal=1:length(Indexes)
    Index=Indexes(animal);

    %% Load metadata and create ouput folder/////////////////////////////////////////////////////////////////////////////
    % Load timing vector
    Timing=csvread([Path filesep 'Chimera' num2str(Index) filesep 'Timing.csv']);

    % Load vector containing the first and last frames tracked
    FramesTracking=csvread([Path filesep 'Chimera' num2str(Index) filesep 'FramesTracking.csv']);

    % Load pixel size (µm)
    Pixelsize=csvread([Path filesep 'Chimera' num2str(Index) filesep 'Pixelsize.csv']);

    % Load registration data
    Registration=csvread([Path filesep 'Chimera' num2str(Index) filesep 'Registration.csv']);

    % Create output folder
    PathOut=[Path filesep 'Chimera' num2str(Index) filesep 'values'];
    mkdir(PathOut);


    %% Determine the range of frames of interest //////////////////////////////////////////////////////////////////////
    % Find the frame corresponding to 2hALS
    temp=abs(Timing-2);
    if (min(temp)<0.5)
        Frame2h=find(temp==min(temp));
    else
        Frame2h=nan;
    end
    % Find the frame corresponding to 12hALS
    temp=abs(Timing-12);
    if (min(temp)<0.5)
        Frame12h=find(temp==min(temp));
    else
        Frame12h=nan;
    end
    % Take the intersect between 2-12h and the tracked period
    if isnan(Frame2h)
        FrameBegin=FramesTracking(1);
    else
        FrameBegin=Frame2h;
    end
    if isnan(Frame12h)
        FrameEnd=FramesTracking(end);
    else
        FrameEnd=Frame12h;
    end

    %% Load manual tracking of hypoblast islands /////////////////////////////////////////////////////////////////////////////
    opts=delimitedTextImportOptions("NumVariables", 7);
    opts.DataLines=[2, Inf];
    opts.Delimiter=",";
    opts.VariableNames=["Trackn", "Slicen", "X", "Y", "Distance", "Velocity", "PixelValue"];
    opts.VariableTypes=["double", "double", "double", "double", "double", "double", "double"];
    opts.ExtraColumnsRule="ignore";
    opts.EmptyLineRule = "read";
    Data=readtable([Path filesep 'Chimera' num2str(Index) filesep 'Tracking.csv'], opts);
    % Extraction of the data from the table
    Trackn=table2array(Data(:,1));
    Slicen=table2array(Data(:,2));
    Xs=table2array(Data(:,3));
    Ys=table2array(Data(:,4));
    % Correct Slicen for the mismatch between manually tracked frames and
    % PIV-tracked frames
    Slicen=Slicen+FramesTracking(1)-1;

    %% Load epiblast PIV data and express values in a Lagrangian fashion /////////////////////////////////////////////////////
    % Initialize the PIV-related vectors/matrices
    Spacing=h5read([Path filesep 'Chimera' num2str(Index) filesep 'data' filesep 'Host-flows.h5'],['/' num2str(1) '/spacing']);
    xmin=h5read([Path filesep 'Chimera' num2str(Index) filesep 'data' filesep 'Host-flows.h5'],['/' num2str(1) '/xmin']);
    xmax=h5read([Path filesep 'Chimera' num2str(Index) filesep 'data' filesep 'Host-flows.h5'],['/' num2str(1) '/xmax']);
    ymin=h5read([Path filesep 'Chimera' num2str(Index) filesep 'data' filesep 'Host-flows.h5'],['/' num2str(1) '/ymin']);
    ymax=h5read([Path filesep 'Chimera' num2str(Index) filesep 'data' filesep 'Host-flows.h5'],['/' num2str(1) '/ymax']);
    BinX=Spacing:Spacing:xmax;
    BinY=Spacing:Spacing:ymax;
    FrameMax=length(FrameBegin:FrameEnd-1);
    SpeedX=nan(length(BinY),length(BinX),FrameMax);
    SpeedY=nan(length(BinY),length(BinX),FrameMax);
    % Fill up the matrices frame after frame
    for frame=FrameBegin:FrameEnd-1
        t=frame-FrameBegin+1;
        % Load h5 file
        dx=h5read([Path filesep 'Chimera' num2str(Index) filesep 'data' filesep 'Host-flows.h5'],['/' num2str(frame) '/dx']);
        dy=h5read([Path filesep 'Chimera' num2str(Index) filesep 'data' filesep 'Host-flows.h5'],['/' num2str(frame) '/dy']);
        x=h5read([Path filesep 'Chimera' num2str(Index) filesep 'data' filesep 'Host-flows.h5'],['/' num2str(frame) '/x']);
        y=h5read([Path filesep 'Chimera' num2str(Index) filesep 'data' filesep 'Host-flows.h5'],['/' num2str(frame) '/y']);

        % Put the values in the appropiate place in the matrix
        for arrow=1:length(x)
            indexX=find(BinX==x(arrow));
            indexY=find(BinY==y(arrow));
            if ~isempty(indexX) && ~isempty(indexY)
                SpeedX(indexY,indexX,t)=dx(arrow);
                SpeedY(indexY,indexX,t)=dy(arrow);
            end
        end
    end

    % Transform epiblast PIV matrices into lagrangian trackings //////////
    % Detect NaN values and interpolate them temporally (in order not to
    % loose trackings at the edge)
    MaskNaNs=isnan(SpeedX);
    VX=inpaintn(SpeedX);
    VY=inpaintn(SpeedY);
    % Initialization of lagrangian position matrices (positions in pixels)
    LagrangePosX=nan(size(VX,2),size(VX,1),size(VX,3)+1);
    LagrangePosY=nan(size(VY,2),size(VY,1),size(VY,3)+1);
    % Timepoint 1 is filled with initial tissue positions
    LagrangePosX(:,:,1)=repmat(BinX,size(VX,2),1);
    LagrangePosY(:,:,1)=repmat(BinY',1,size(VX,1));

    % Filling lagrangian position matrices in time
    tic
    for t=1:size(VX,3)
        % Interpolate the PIV grid to the size of the image to get a speed value for each pixel of the image
        temp=VX(:,:,t);
        [X,Y]=meshgrid(BinY,BinX);
        [Xq,Yq]=meshgrid(1:ymax,1:xmax);
        InterpVX=interp2(X,Y,temp,Xq,Yq);
        temp=VY(:,:,t);
        InterpVY=interp2(X,Y,temp,Xq,Yq);

        % Pixel by pixel filling up
        for x=1:size(VX,2)
            for y=1:size(VX,1)
                % Index for sampling VX and VY to calculate
                % displacement (in pixels)
                xrep=round(LagrangePosX(y,x,t));
                yrep=round(LagrangePosY(y,x,t));
                % If index are notvalid (NaN value, outside the field),
                % fill the lagrangian position matrice with a NaN
                if isnan(xrep) || isnan(yrep) || xrep>xmax || yrep>ymax || xrep<1 || yrep<1
                    LagrangePosX(y,x,t+1)=nan;
                    LagrangePosY(y,x,t+1)=nan;
                else
                    % If the index is valid, update the position by adding the recorded displacement (in pixels)
                    LagrangePosX(y,x,t+1)=LagrangePosX(y,x,t)+InterpVX(yrep,xrep);
                    LagrangePosY(y,x,t+1)=LagrangePosY(y,x,t)+InterpVY(yrep,xrep);
                end
            end
        end
    end
    toc

    % Put back NaN values
    Mask=MaskNaNs(:,:,1);
    for t=1:size(LagrangePosX,3)
        temp=LagrangePosX(:,:,t);
        temp(Mask)=NaN;
        LagrangePosX(:,:,t)=temp;
        temp=LagrangePosY(:,:,t);
        temp(Mask)=NaN;
        LagrangePosY(:,:,t)=temp;
    end

    % Interpolate so as to have a value for each pixel
    [X,Y,T]=meshgrid(BinY,BinX,1:size(LagrangePosX,3));
    [Xq,Yq,Tq]=meshgrid(1:ymax,1:xmax,1:size(LagrangePosX,3));
    InterpLagrangePosX=interp3(X,Y,T,LagrangePosX,Xq,Yq,Tq);
    InterpLagrangePosY=interp3(X,Y,T,LagrangePosY,Xq,Yq,Tq);
    [X,Y,T]=meshgrid(BinY,BinX,1:size(VX,3));
    [Xq,Yq,Tq]=meshgrid(1:ymax,1:xmax,1:size(VX,3));
    InterpVX=interp3(X,Y,T,VX,Xq,Yq,Tq);
    InterpVY=interp3(X,Y,T,VY,Xq,Yq,Tq);   

    %% Calculate differential motion (ref timing 4h) ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % Identify total number of tracks
    NumberOfTracks=max(Trackn);
    % Identify timestep to define smoothing radius
    Timestep=Timing(2)-Timing(1);
    SmoothingRadius=round(SmoothingRadiusT/Timestep);
    % Identify the frame corresponding to 4h
    temp=abs(Timing-4);
    if (min(temp)<0.5)
        Frame4h=find(temp==min(temp));
    else
        Frame4h=nan;
    end
    % Initialization of vectors for saving differential motion and hypo/epi
    % motions
    TracknDiff=[];
    SlicenDiff=[];
    XsDiff=[];
    YsDiff=[];
    XsSmoothed=[];
    YsSmoothed=[];
    XsEpi=[];
    YsEpi=[];
    XsEpi4h=[];
    YsEpi4h=[];

    % For each tracked hypoblast island
    for track=1:NumberOfTracks
        % Isolate track
        index=find(Trackn==track);
        XTrack=Xs(index);
        YTrack=Ys(index);
        SliceTrack=Slicen(index);
        % Smooth the track
        XTrack=smooth(XTrack,SmoothingRadius);
        YTrack=smooth(YTrack,SmoothingRadius);
        % Count the number of points of interest to trim all the vectors
        index=find(SliceTrack>=FrameBegin & SliceTrack<=FrameEnd);
        % Trim the vectors to cover only the timings of interest
        XTrack=XTrack(index);
        YTrack=YTrack(index);
        SliceTrack=SliceTrack(index);
        % Initialize vectors storing differential motion
        DiffXTrack=nan(1,length(SliceTrack));
        DiffYTrack=nan(1,length(SliceTrack));
        DiffXTrack(1)=XTrack(1);
        DiffYTrack(1)=YTrack(1);
        % Initialize vectors storing epiblast motion
        EpiXTrack=nan(1,length(SliceTrack));
        EpiYTrack=nan(1,length(SliceTrack));
        EpiXTrack(1)=XTrack(1);
        EpiYTrack(1)=YTrack(1);
        % Initialize vectors storing epiblast motion relative to 4h
        % position
        Slice4h=find(SliceTrack==Frame4h);
        EpiXTrack4h=nan(1,length(SliceTrack));
        EpiYTrack4h=nan(1,length(SliceTrack));
        EpiXTrack4h(Slice4h)=XTrack(Slice4h);
        EpiYTrack4h(Slice4h)=YTrack(Slice4h);
        % For each timepoint of interest
        for slice=1:length(index)-1
            % Index for matrices selection
            t=SliceTrack(slice)-FrameBegin+1;
            % Extract the motion of the hypoblast island from t to (t+1)
            dxhypo=XTrack(slice+1)-XTrack(slice);
            dyhypo=YTrack(slice+1)-YTrack(slice);
            % Extract the motion of the epiblast from t to (t+1) at the
            % same location of the hypoblast at time t
            if round(XTrack(slice))<=xmax && round(YTrack(slice))<=ymax &&  round(XTrack(slice))>=1 && round(YTrack(slice))>=1 && ~isnan(XTrack(slice)) && ~isnan(YTrack(slice))
                dxepi=InterpVX(round(YTrack(slice)),round(XTrack(slice)),t);
                dyepi=InterpVY(round(YTrack(slice)),round(XTrack(slice)),t);
                % Update the track with the motion of the hypoblast minus
                % the motion of the epiblast
                DiffXTrack(slice+1)=dxhypo-dxepi+DiffXTrack(slice);
                DiffYTrack(slice+1)=dyhypo-dyepi+DiffYTrack(slice);
            else
                DiffXTrack(slice+1)=nan;
                DiffYTrack(slice+1)=nan;
            end
            % Extract motion of the epiblast relative to the first time
            % point
            if round(EpiXTrack(slice))<=xmax && round(EpiYTrack(slice))<=ymax &&  round(EpiXTrack(slice))>=1 && round(EpiYTrack(slice))>=1 && ~isnan(EpiXTrack(slice)) && ~isnan(EpiYTrack(slice))
                dxepi=InterpLagrangePosX(round(EpiYTrack(1)),round(EpiXTrack(1)),t+1)-InterpLagrangePosX(round(EpiYTrack(1)),round(EpiXTrack(1)),t);
                dyepi=InterpLagrangePosY(round(EpiYTrack(1)),round(EpiXTrack(1)),t+1)-InterpLagrangePosY(round(EpiYTrack(1)),round(EpiXTrack(1)),t);
                EpiXTrack(slice+1)=dxepi+EpiXTrack(slice);
                EpiYTrack(slice+1)=dyepi+EpiYTrack(slice);
            else
                EpiXTrack(slice+1)=nan;
                EpiYTrack(slice+1)=nan;
            end
            % Extract motion of the epiblast relative to position at 4h
            if ~isempty(Slice4h) && slice>=Slice4h && round(EpiXTrack4h(slice))<=xmax && round(EpiYTrack4h(slice))<=ymax &&  round(EpiXTrack4h(slice))>=1 && round(EpiYTrack4h(slice))>=1 && ~isnan(EpiXTrack4h(slice)) && ~isnan(EpiYTrack4h(slice))
                dxepi=InterpLagrangePosX(round(EpiYTrack4h(Slice4h)),round(EpiXTrack4h(Slice4h)),t+1)-InterpLagrangePosX(round(EpiYTrack4h(Slice4h)),round(EpiXTrack4h(Slice4h)),t);
                dyepi=InterpLagrangePosY(round(EpiYTrack4h(Slice4h)),round(EpiXTrack4h(Slice4h)),t+1)-InterpLagrangePosY(round(EpiYTrack4h(Slice4h)),round(EpiXTrack4h(Slice4h)),t);
                EpiXTrack4h(slice+1)=dxepi+EpiXTrack4h(slice);
                EpiYTrack4h(slice+1)=dyepi+EpiYTrack4h(slice);
            elseif ~isempty(Slice4h) && slice<=Slice4h
            else
                EpiXTrack4h(slice+1)=nan;
                EpiYTrack4h(slice+1)=nan;
            end
        end
        % Store the track for saving coordinates
        TracknDiff=cat(2,TracknDiff,track*ones(1,length(SliceTrack)));
        SlicenDiff=cat(2,SlicenDiff,SliceTrack');
        XsDiff=cat(2,XsDiff,DiffXTrack);
        YsDiff=cat(2,YsDiff,DiffYTrack);
        XsSmoothed=cat(2,XsSmoothed,XTrack');
        YsSmoothed=cat(2,YsSmoothed,YTrack');
        XsEpi=cat(2,XsEpi,EpiXTrack);
        YsEpi=cat(2,YsEpi,EpiYTrack); 
        XsEpi4h=cat(2,XsEpi4h,EpiXTrack4h);
        YsEpi4h=cat(2,YsEpi4h,EpiYTrack4h);
    end
    
    % Save differential motion vectors for Fiji plotting
    writematrix(TracknDiff,[PathOut filesep 'Trackn.txt']);
    writematrix(SlicenDiff,[PathOut filesep 'Slicen.txt']);
    writematrix(XsDiff,[PathOut filesep 'XsDiff.txt']);
    writematrix(YsDiff,[PathOut filesep 'YsDiff.txt']);
    writematrix(XsSmoothed,[PathOut filesep 'XsSmoothed.txt']);
    writematrix(YsSmoothed,[PathOut filesep 'YsSmoothed.txt']);
    writematrix(XsEpi,[PathOut filesep 'XsEpi.txt']);
    writematrix(YsEpi,[PathOut filesep 'YsEpi.txt']);
    % Save useful data for plotting
    writematrix([xmax ymax],[PathOut filesep 'ImageSize.txt']);
    writematrix(max(SlicenDiff),[PathOut filesep 'MaxFrame.txt']);
    writematrix(max(TracknDiff),[PathOut filesep 'NumberTracks.txt']);

 %% Plot average differential motion/ hypoblast motion vector maps /////////////////////////////////////////////////////////////////////////
    % Extract displacements from 4 to 12h
    HypodX4_12h=[];
    HypodY4_12h=[];
    DiffX4_12h=[];
    DiffY4_12h=[];
    for track=1:NumberOfTracks
        index=find(TracknDiff==track);
        XTrack=XsSmoothed(index);
        YTrack=YsSmoothed(index);
        DiffXTrack=XsDiff(index);
        DiffYTrack=YsDiff(index);
        SliceTrack=SlicenDiff(index);
        if ~isempty(find(SliceTrack==Frame4h,1)) && ~isempty(find(SliceTrack==Frame12h,1)) 
            index=find(SliceTrack==Frame4h);
            index2=find(SliceTrack==Frame12h);            
            HypodX4_12h=cat(2,HypodX4_12h,XTrack(index2)-XTrack(index));
            HypodY4_12h=cat(2,HypodY4_12h,YTrack(index2)-YTrack(index));
            DiffX4_12h=cat(2,DiffX4_12h,DiffXTrack(index2)-DiffXTrack(index));
            DiffY4_12h=cat(2,DiffY4_12h,DiffYTrack(index2)-DiffYTrack(index));
        end
    end

    % Epiblast motion
    EpidX4_12h=[];
    EpidY4_12h=[];
    X4_12h=[];
    Y4_12h=[];
    % For each track, fill the matrices
    for track=1:NumberOfTracks
        % Isolate track
        index=find(TracknDiff==track);
        XTrack=XsSmoothed(index);
        YTrack=YsSmoothed(index);
        DiffXTrack=XsDiff(index);
        DiffYTrack=YsDiff(index);
        EpiXTrack=XsEpi4h(index);
        EpiYTrack=YsEpi4h(index);
        SliceTrack=SlicenDiff(index);
        if ~isempty(find(SliceTrack==Frame4h,1)) && ~isempty(find(SliceTrack==Frame12h,1))
            index=find(SliceTrack==Frame4h);
            index3=find(SliceTrack==Frame12h);   
            X4_12h=cat(2,X4_12h,XTrack(index));
            Y4_12h=cat(2,Y4_12h,YTrack(index));
            EpidX4_12h=cat(2,EpidX4_12h,EpiXTrack(index3)-EpiXTrack(index));
            EpidY4_12h=cat(2,EpidY4_12h,EpiYTrack(index3)-EpiYTrack(index));
        end
    end

    %% Registration of the data to generate an archetype map
    % Put the landmark point at the origin
    X4_12h=X4_12h-Registration(1);
    Y4_12h=Y4_12h-Registration(2);

    % Normalize for the diameter of the embryo
    X4_12h=100*X4_12h/Registration(4);
    Y4_12h=100*Y4_12h/Registration(4);

    % Transform position and values into polar coordinates and apply the
    % rotation
    [theta,rho]=cart2pol(X4_12h,Y4_12h);
    [Hypotheta,Hyporho]=cart2pol(HypodX4_12h,HypodY4_12h);
    [Difftheta,Diffrho]=cart2pol(DiffX4_12h,DiffY4_12h);
    [Epitheta,Epirho]=cart2pol(EpidX4_12h,EpidY4_12h);
    theta=theta+deg2rad(Registration(3)+90);
    Hypotheta=Hypotheta+deg2rad(Registration(3)+90);
    Difftheta=Difftheta+deg2rad(Registration(3)+90);
    Epitheta=Epitheta+deg2rad(Registration(3)+90);
    % After rotation, put back in cartesien coordinates
    [X4_12hRotated,Y4_12hRotated]=pol2cart(theta,rho);
    [HypodX4_12hRotated,HypodY4_12hRotated]=pol2cart(Hypotheta,Hyporho);
    [DiffX4_12hRotated,DiffY4_12hRotated]=pol2cart(Difftheta,Diffrho);
    [EpidX4_12hRotated,EpidY4_12hRotated]=pol2cart(Epitheta,Epirho);

    % Normalize the displacement by the height of the animal (% of height)
    HypodX4_12hRotated=HypodX4_12hRotated*100/Registration(4);
    HypodY4_12hRotated=HypodY4_12hRotated*100/Registration(4);
    DiffX4_12hRotated=DiffX4_12hRotated*100/Registration(4);
    DiffY4_12hRotated=DiffY4_12hRotated*100/Registration(4);
    EpidX4_12hRotated=EpidX4_12hRotated*100/Registration(4);
    EpidY4_12hRotated=EpidY4_12hRotated*100/Registration(4);

    % Save values
    PoolX=cat(2,PoolX,X4_12hRotated);
    PoolY=cat(2,PoolY,Y4_12hRotated);
    PoolAnimal=cat(2,PoolAnimal,animal*ones(1,length(Y4_12hRotated)));
    PoolHypoX=cat(2,PoolHypoX,HypodX4_12hRotated);
    PoolHypoY=cat(2,PoolHypoY,HypodY4_12hRotated);
    PoolDiffX=cat(2,PoolDiffX,DiffX4_12hRotated);
    PoolDiffY=cat(2,PoolDiffY,DiffY4_12hRotated);
    PoolEpiX=cat(2,PoolEpiX,EpidX4_12hRotated);
    PoolEpiY=cat(2,PoolEpiY,EpidY4_12hRotated);    

end

% % Saving (this section is commented not to overwrite the full dataset, as
% here only one example of the raw data is provided
% save([PathData filesep 'PoolX'],'PoolX');
% save([PathData filesep 'PoolY'],'PoolY');
% save([PathData filesep 'PoolAnimal'],'PoolAnimal');
% save([PathData filesep 'PoolHypoX'],'PoolHypoX');
% save([PathData filesep 'PoolHypoY'],'PoolHypoY');
% save([PathData filesep 'PoolEpiX'],'PoolEpiX');
% save([PathData filesep 'PoolEpiY'],'PoolEpiY');
% save([PathData filesep 'PoolDiffX'],'PoolDiffX');
% save([PathData filesep 'PoolDiffY'],'PoolDiffY');

%% Plotting

% Load data
load([PathData filesep 'PoolX.mat']);
load([PathData filesep 'PoolY.mat']);
load([PathData filesep 'PoolAnimal.mat']);
load([PathData filesep 'PoolHypoX.mat']);
load([PathData filesep 'PoolHypoY.mat']);
load([PathData filesep 'PoolEpiX.mat']);
load([PathData filesep 'PoolEpiY.mat']);
load([PathData filesep 'PoolDiffX.mat']);
load([PathData filesep 'PoolDiffY.mat']);

% Plot all tracks (pooling all the registered movies)
factor=0.15;
figure(1)
quiver(PoolX,PoolY,PoolEpiX*factor,PoolEpiY*factor,'AutoScale','off','ShowArrowHead','on','Color','green')
daspect([1 1 1])
title('All epiblast tracks')
figure(2)
quiver(PoolX,PoolY,PoolHypoX*factor,PoolHypoY*factor,'AutoScale','off','ShowArrowHead','on','Color','magenta')
daspect([1 1 1])
title('All hypoblast tracks')
figure(3)
quiver(PoolX,PoolY,PoolDiffX*factor,PoolDiffY*factor,'AutoScale','off','ShowArrowHead','on','Color','blue')
daspect([1 1 1])
title('All differential motion tracks')

% Spatial averaging
step=10;
Xbin=-80:step:80;
Ybin=-60:step:140;
AverageDiffX=nan(length(Xbin),length(Ybin));
AverageDiffY=nan(length(Xbin),length(Ybin));
AverageHypoX=nan(length(Xbin),length(Ybin));
AverageHypoY=nan(length(Xbin),length(Ybin));
AverageEpiX=nan(length(Xbin),length(Ybin));
AverageEpiY=nan(length(Xbin),length(Ybin));
for x=1:length(Xbin)
    indexX=find(PoolX>=Xbin(x)-step/2 & PoolX<=Xbin(x)+step/2);
    for y=1:length(Ybin)
        indexY=find(PoolY>=Ybin(y)-step/2 & PoolY<=Ybin(y)+step/2);
        index=intersect(indexX,indexY);
        if length(index)>=2
            AverageDiffX(x,y)=mean(PoolDiffX(index));
            AverageDiffY(x,y)=mean(PoolDiffY(index));
            AverageHypoX(x,y)=mean(PoolHypoX(index));
            AverageHypoY(x,y)=mean(PoolHypoY(index));
            AverageEpiX(x,y)=mean(PoolEpiX(index));
            AverageEpiY(x,y)=mean(PoolEpiY(index));
        end
    end
end

% Plot average velocity fields
factor=0.5;
figure(4)
quiver(repmat(Xbin',1,length(Ybin)),repmat(Ybin,length(Xbin),1),AverageEpiX*factor,AverageEpiY*factor,'AutoScale','off','ShowArrowHead','on','Color','green','LineWidth',1.5)
daspect([1 1 1])
title('Epiblast velocity field')
figure(5)
quiver(repmat(Xbin',1,length(Ybin)),repmat(Ybin,length(Xbin),1),AverageHypoX*factor,AverageHypoY*factor,'AutoScale','off','ShowArrowHead','on','Color','magenta','LineWidth',1.5)
daspect([1 1 1])
title('Hypoblast velocity field')
figure(6)
quiver(repmat(Xbin',1,length(Ybin)),repmat(Ybin,length(Xbin),1),AverageDiffX*factor,AverageDiffY*factor,'AutoScale','off','ShowArrowHead','on','LineWidth',1.5)
daspect([1 1 1])
title('Differential motion velocity field')

% Grid Epiblast
EpiX=inpaintn(AverageEpiX);
EpiY=inpaintn(AverageEpiY);
% Smoothing
EpiX=smoothdata2(EpiX,"movmean",3);
EpiY=smoothdata2(EpiY,"movmean",3);
% Lines drawing
figure(7)
hold on
for y=1:size(EpiX,2)
    X=Xbin+squeeze(EpiX(:,y))';
    Y=Ybin(y)+squeeze(EpiY(:,y))';
    Vect=Sampling(:,y);
    index=find(Vect==1);
    plot(X(index),Y(index),'Color','green');
end
for x=1:size(EpiX,1)
    X=Xbin(x)+squeeze(EpiX(x,:))';
    Y=Ybin'+squeeze(EpiY(x,:))';
    Vect=Sampling(x,:);
    index=find(Vect==1);
    plot(X(index),Y(index),'Color','green');
end
daspect([1 1 1])
hold off
title('Grid epiblast')

% Grid Hypoblast
HypoX=inpaintn(AverageHypoX);
HypoY=inpaintn(AverageHypoY);
% Smoothing
HypoX=smoothdata2(HypoX,"movmean",3);
HypoY=smoothdata2(HypoY,"movmean",3);
% Lines drawing
figure(8)
hold on
for y=1:size(HypoX,2)
    X=Xbin+squeeze(HypoX(:,y))';
    Y=Ybin(y)+squeeze(HypoY(:,y))';
    Vect=Sampling(:,y);
    index=find(Vect==1);
    plot(X(index),Y(index),'Color','magenta');
end
for x=1:size(DiffX,1)
    X=Xbin(x)+squeeze(HypoX(x,:))';
    Y=Ybin'+squeeze(HypoY(x,:))';
    Vect=Sampling(x,:);
    index=find(Vect==1);
    plot(X(index),Y(index),'Color','magenta');
end
daspect([1 1 1])
hold off
title('Grid hypoblast')


% Grid Difference
Sampling=double(~isnan(AverageDiffX));
SamplingSmoothed=smoothdata2(Sampling,"movmedian",3);
Sampling=double(SamplingSmoothed==1);
% Interpolation
DiffX=inpaintn(AverageDiffX);
DiffY=inpaintn(AverageDiffY);
% Smoothing
DiffX=smoothdata2(DiffX,"movmean",3);
DiffY=smoothdata2(DiffY,"movmean",3);
% Putting back nans
DiffX(Sampling==0)=nan;
DiffY(Sampling==0)=nan;
% Lines drawing
figure(9)
hold on
for y=1:size(DiffX,2)
    X=Xbin+squeeze(DiffX(:,y))';
    Y=Ybin(y)+squeeze(DiffY(:,y))';
    Vect=Sampling(:,y);
    index=find(Vect==1);
    plot(X(index),Y(index),'Color','blue');
end
for x=1:size(DiffX,1)
    X=Xbin(x)+squeeze(DiffX(x,:))';
    Y=Ybin'+squeeze(DiffY(x,:))';
    Vect=Sampling(x,:);
    index=find(Vect==1);
    plot(X(index),Y(index),'Color','blue');
end
daspect([1 1 1])
hold off
xlim([-100 100])
ylim([-60 120])
title('Deformation grid differential motion')