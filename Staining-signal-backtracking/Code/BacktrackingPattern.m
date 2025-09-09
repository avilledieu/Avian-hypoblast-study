%% BacktrackingPattern
% Bactrack a pattern of interest (aligned with the last frame of a movie)
% to 4h

clearvars
close all

%% Parameters
% Path where to find embryo subfolders
Path='\\gaia.pasteur.fr\MVS_DATA1\Aurelien\9_Fixed-samples\24-06-06_SUMMARY-Nodal-FoxA2-patterns-in-time';
% Timings to consider
timing='8h';
% Indexes of the embryos associated with each timing
Index=1:7;

%% Code

% Load average residual hypoblast motion after substraction of epiblast
% motion from 4h to 8h (calculated using 2-color chimeras)
load([Path filesep 'DiffEpiHypo.mat']);

for embryo=Index
    tic
    %% Calculation using PIV data to determine bactracking coordinates

    % Read the PIV metadata (epiblast)
    Spacing=h5read([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'data' filesep 'Max-flows.h5'],['/' num2str(1) '/spacing']);
    xmin=h5read([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'data' filesep 'Max-flows.h5'],['/' num2str(1) '/xmin']);
    xmax=h5read([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'data' filesep 'Max-flows.h5'],['/' num2str(1) '/xmax']);
    ymin=h5read([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'data' filesep 'Max-flows.h5'],['/' num2str(1) '/ymin']);
    ymax=h5read([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'data' filesep 'Max-flows.h5'],['/' num2str(1) '/ymax']);
    BinX=Spacing:Spacing:xmax;
    BinY=Spacing:Spacing:ymax;

    % Read spatial and temporal alignment data
    Timing=csvread([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'Timing.csv']);
    CoordinatesCentersRotation=csvread([Path filesep timing '_' num2str(embryo) filesep 'CoordinatesCentersRotation.csv']);

    % Determine 4h frame
    temp=abs(Timing-4);
    if (min(temp)<0.5)
        Frame4h=find(temp==min(temp));
    else
        Frame4h=nan;
    end
    FrameMax=length(Timing);

    % Initialization of the PIV-related vectors/matrices
    SpeedX=nan(length(BinY),length(BinX),FrameMax-Frame4h-1);
    SpeedY=nan(length(BinY),length(BinX),FrameMax-Frame4h-1);

    % Fill up the PIV matrices frame after frame (from 4h to 8h, epiblast)
    for frame=Frame4h:FrameMax-2
        t=frame+1-Frame4h;
        % Load h5 file
        dx=h5read([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'data' filesep 'Max-flows.h5'],['/' num2str(frame) '/dx']);
        dy=h5read([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'data' filesep 'Max-flows.h5'],['/' num2str(frame) '/dy']);
        x=h5read([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'data' filesep 'Max-flows.h5'],['/' num2str(frame) '/x']);
        y=h5read([Path filesep timing '_' num2str(embryo) filesep 'Live' filesep 'data' filesep 'Max-flows.h5'],['/' num2str(frame) '/y']);

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

    
    % Reshape the 'hypoblast minus epiblast motion' field to make it correspond with the
    % last time point of hypoblast movie (make tilt and center of rotations
    % correspond)
    % The centers were calculated on a 600X600pixels degraded image. First,
    % we put them back in original image value (xmax,ymax)
    CenterX1=CoordinatesCentersRotation(1)/600*xmax;
    CenterY1=CoordinatesCentersRotation(2)/600*ymax;
    CenterX2=CoordinatesCentersRotation(3)/600*xmax;
    CenterY2=CoordinatesCentersRotation(4)/600*ymax;
    % Calculate tilt of the epiblast image
    tilt=atan((CenterY2-CenterY1)/(CenterX2-CenterX1));
    % Calculate the midpoint between the two center in the image
    MidpointX=(CenterX1+CenterX2)/2;
    MidpointY=(CenterY1+CenterY2)/2;
    % Calculate the distance between the two centers of rotation (in
    % pixels)
    EmbryoLength=sqrt((CenterY2-CenterY1)^2+(CenterX2-CenterX1)^2);
    % Make the midpoint of the archetype coordinates match the midpoint of
    % the image (midpoint of the archetype coordinates = 0, 35)
    Xbintemp=DiffEpiHypo.Xbin;
    Ybintemp=DiffEpiHypo.Ybin;
    Xbin=repmat(Xbintemp',1,length(Ybintemp));
    Ybin=repmat(Ybintemp,length(Xbintemp),1);
    Xbin=Xbin;
    Ybin=Ybin-35;
    % Make pixel sizes match in coordinates and displacement values (in the
    % archetype, embryo length = 70
    Xbin=Xbin/70*EmbryoLength;
    Ybin=Ybin/70*EmbryoLength;
    SpeedDiffX=DiffEpiHypo.AverageDiffX4_8h/70*EmbryoLength;
    SpeedDiffY=DiffEpiHypo.AverageDiffY4_8h/70*EmbryoLength;
    Xbin=Xbin+MidpointX;
    Ybin=Ybin+MidpointY;
    % Apply to the archetype coordinates a rotation of center Midpoint and
    % of angle -tilt
    XbinRotat=(Xbin-MidpointX)*cos(-tilt)-(Ybin-MidpointY)*sin(-tilt)+MidpointX;
    YbinRotat=(Xbin-MidpointX)*sin(-tilt)+(Ybin-MidpointY)*cos(-tilt)+MidpointY;
    % Apply a rotation of tilt to the displacements
    SpeedDiffXRotat=SpeedDiffX*cos(-tilt)-SpeedDiffY*sin(-tilt);
    SpeedDiffYRotat=SpeedDiffX*sin(-tilt)+SpeedDiffY*cos(-tilt);
    % % Visual check of the rotation
    % figure()
    % subplot(1,2,1)
    % quiver(repmat(Xbintemp',1,length(Ybintemp)),repmat(Ybintemp,length(Xbintemp),1),SpeedDiffX,SpeedDiffY,'AutoScale','on','ShowArrowHead','on','Color','blue','LineWidth',1.5)
    % title('Before rotation and rescaling')
    % daspect([1 1 1])
    % subplot(1,2,2)
    % quiver(XbinRotat,YbinRotat,SpeedDiffXRotat,SpeedDiffYRotat,'AutoScale','on','ShowArrowHead','on','Color','green','LineWidth',1.5)
    % title('After rotation and rescaling')
    % daspect([1 1 1])

    % Interpolate the hypoblast flow field so that it matchs PIV boxes
    X=YbinRotat;
    Y=XbinRotat;
    [Xq,Yq]=meshgrid(BinY,BinX);
    SpeedDiffXInterp=griddata(X,Y,SpeedDiffXRotat,Xq,Yq);
    SpeedDiffYInterp=griddata(X,Y,SpeedDiffYRotat,Xq,Yq);
    SpeedDiffXInpaint=inpaintn(SpeedDiffXInterp);
    SpeedDiffYInpaint=inpaintn(SpeedDiffYInterp);
    % quiver(Xq,Yq,SpeedDiffXInpaint,SpeedDiffYInpaint,'AutoScale','on','ShowArrowHead','on','Color','green','LineWidth',1.5)

    % Add the 'hypoblast minus epiblast' flows as an extra
    % timepoint for hypoblast backtracking
    SpeedXHypo=cat(3,SpeedX,SpeedDiffXInpaint);
    SpeedYHypo=cat(3,SpeedY,SpeedDiffYInpaint);

    %% Transform PIV matrices into lagrangian trackings (epiblast) ////////
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

    % Extract lagrangian position of last frame of the movie
    LagrangePosX=InterpLagrangePosX(:,:,end);
    LagrangePosY=InterpLagrangePosY(:,:,end);
    LagrangePosX=round(LagrangePosX);
    LagrangePosY=round(LagrangePosY);

 %% Transform PIV matrices into lagrangian trackings (hypoblast) /////////
    % Detect NaN values and interpolate them temporally (in order not to
    % loose trackings at the edge)
    MaskNaNs=isnan(SpeedXHypo);
    VX=inpaintn(SpeedXHypo);
    VY=inpaintn(SpeedYHypo);
    % Initialization of lagrangian position matrices (positions in pixels)
    LagrangePosXHypo=nan(size(VX,2),size(VX,1),size(VX,3)+1);
    LagrangePosYHypo=nan(size(VY,2),size(VY,1),size(VY,3)+1);
    % Timepoint 1 is filled with initial tissue positions
    LagrangePosXHypo(:,:,1)=repmat(BinX,size(VX,2),1);
    LagrangePosYHypo(:,:,1)=repmat(BinY',1,size(VX,1));

    % Filling lagrangian position matrices in time
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
                xrep=round(LagrangePosXHypo(y,x,t));
                yrep=round(LagrangePosYHypo(y,x,t));
                % If index are notvalid (NaN value, outside the field),
                % fill the lagrangian position matrice with a NaN
                if isnan(xrep) || isnan(yrep) || xrep>xmax || yrep>ymax || xrep<1 || yrep<1
                    LagrangePosXHypo(y,x,t+1)=nan;
                    LagrangePosYHypo(y,x,t+1)=nan;
                else
                    % If the index is valid, update the position by adding the recorded displacement (in pixels)
                    LagrangePosXHypo(y,x,t+1)=LagrangePosXHypo(y,x,t)+InterpVX(yrep,xrep);
                    LagrangePosYHypo(y,x,t+1)=LagrangePosYHypo(y,x,t)+InterpVY(yrep,xrep);
                end
            end
        end
    end

    % Put back NaN values
    Mask=MaskNaNs(:,:,1);
    for t=1:size(LagrangePosX,3)
        temp=LagrangePosXHypo(:,:,t);
        temp(Mask)=NaN;
        LagrangePosXHypo(:,:,t)=temp;
        temp=LagrangePosYHypo(:,:,t);
        temp(Mask)=NaN;
        LagrangePosYHypo(:,:,t)=temp;
    end

    % Interpolate so as to have a value for each pixel
    [X,Y,T]=meshgrid(BinY,BinX,1:size(LagrangePosXHypo,3));
    [Xq,Yq,Tq]=meshgrid(1:ymax,1:xmax,1:size(LagrangePosXHypo,3));
    InterpLagrangePosXHypo=interp3(X,Y,T,LagrangePosXHypo,Xq,Yq,Tq);
    InterpLagrangePosYHypo=interp3(X,Y,T,LagrangePosYHypo,Xq,Yq,Tq);

    % Extract lagrangian position of last frame of the movie
    LagrangePosXHypo=InterpLagrangePosXHypo(:,:,end);
    LagrangePosYHypo=InterpLagrangePosYHypo(:,:,end);
    LagrangePosXHypo=round(LagrangePosXHypo);
    LagrangePosYHypo=round(LagrangePosYHypo);

    %% Apply backtracking to the patterns

    % Open NODAL_dorsal
    Image=imread([Path filesep timing '_' num2str(embryo) filesep 'NODAL_dorsal.tif']);
    % Filling up the backtracked image
    Backtracked=nan(size(Image));
    for x=1:size(Image,2)
        for y=1:size(Image,1)
            X=LagrangePosX(y,x);
            Y=LagrangePosY(y,x);
            if X>=1 && X<=size(Image,2) && Y>=1 && Y<=size(Image,1) && ~isnan(X) && ~isnan(Y)
                Backtracked(y,x)=Image(Y,X);
            end
        end
    end
    % Save backtracked image
    imwrite(uint16(Backtracked),[Path filesep timing '_' num2str(embryo) filesep 'NODAL_dorsal(backtracked4h).tif']);

    % Open NODAL_ventral
    Image=imread([Path filesep timing '_' num2str(embryo) filesep 'NODAL_ventral.tif']);
    % Filling up the backtracked image
    Backtracked=nan(size(Image));
    for x=1:size(Image,2)
        for y=1:size(Image,1)
            X=LagrangePosXHypo(y,x);
            Y=LagrangePosYHypo(y,x);
            if X>=1 && X<=size(Image,2) && Y>=1 && Y<=size(Image,1) && ~isnan(X) && ~isnan(Y)
                Backtracked(y,x)=Image(Y,X);
            end
        end
    end
    % Save backtracked image
    imwrite(uint16(Backtracked),[Path filesep timing '_' num2str(embryo) filesep 'NODAL_ventral(backtracked4h).tif']);

    % Open binarized NODAL_dorsal
    Image=imread([Path filesep timing '_' num2str(embryo) filesep 'NODAL_dorsal(binarized).tif']);
    % Filling up the backtracked image
    Backtracked=nan(size(Image));
    for x=1:size(Image,2)
        for y=1:size(Image,1)
            X=LagrangePosX(y,x);
            Y=LagrangePosY(y,x);
            if X>=1 && X<=size(Image,2) && Y>=1 && Y<=size(Image,1) && ~isnan(X) && ~isnan(Y)
                Backtracked(y,x)=Image(Y,X);
            end
        end
    end
    % Save backtracked image
    imwrite(uint16(Backtracked),[Path filesep timing '_' num2str(embryo) filesep 'NODAL_dorsal(binarized-backtracked4h).tif']);

    % Open binarized NODAL_ventral
    Image=imread([Path filesep timing '_' num2str(embryo) filesep 'NODAL_ventral(binarized).tif']);
    % Filling up the backtracked image
    Backtracked=nan(size(Image));
    for x=1:size(Image,2)
        for y=1:size(Image,1)
            X=LagrangePosX(y,x);
            Y=LagrangePosY(y,x);
            if X>=1 && X<=size(Image,2) && Y>=1 && Y<=size(Image,1) && ~isnan(X) && ~isnan(Y)
                Backtracked(y,x)=Image(Y,X);
            end
        end
    end
    % Save backtracked image
    imwrite(uint16(Backtracked),[Path filesep timing '_' num2str(embryo) filesep 'NODAL_ventral(binarized-backtracked4h).tif']);



    % Open FoxA2_dorsal
    Image=imread([Path filesep timing '_' num2str(embryo) filesep 'FoxA2_dorsal.tif']);
    % Filling up the backtracked image
    Backtracked=nan(size(Image));
    for x=1:size(Image,2)
        for y=1:size(Image,1)
            X=LagrangePosX(y,x);
            Y=LagrangePosY(y,x);
            if X>=1 && X<=size(Image,2) && Y>=1 && Y<=size(Image,1) && ~isnan(X) && ~isnan(Y)
                Backtracked(y,x)=Image(Y,X);
            end
        end
    end
    % Save backtracked image
    imwrite(uint16(Backtracked),[Path filesep timing '_' num2str(embryo) filesep 'FoxA2_dorsal(backtracked4h).tif']);

    % Open FoxA2_ventral
    Image=imread([Path filesep timing '_' num2str(embryo) filesep 'FoxA2_ventral.tif']);
    % Filling up the backtracked image
    Backtracked=nan(size(Image));
    for x=1:size(Image,2)
        for y=1:size(Image,1)
            X=LagrangePosXHypo(y,x);
            Y=LagrangePosYHypo(y,x);
            if X>=1 && X<=size(Image,2) && Y>=1 && Y<=size(Image,1) && ~isnan(X) && ~isnan(Y)
                Backtracked(y,x)=Image(Y,X);
            end
        end
    end
    % Save backtracked image
    imwrite(uint16(Backtracked),[Path filesep timing '_' num2str(embryo) filesep 'FoxA2_ventral(backtracked4h).tif']);

    % Open binarized FoxA2_ventral
    Image=imread([Path filesep timing '_' num2str(embryo) filesep 'FoxA2_ventral(binarized).tif']);
    % Filling up the backtracked image
    Backtracked=nan(size(Image));
    for x=1:size(Image,2)
        for y=1:size(Image,1)
            X=LagrangePosXHypo(y,x);
            Y=LagrangePosYHypo(y,x);
            if X>=1 && X<=size(Image,2) && Y>=1 && Y<=size(Image,1) && ~isnan(X) && ~isnan(Y)
                Backtracked(y,x)=Image(Y,X);
            end
        end
    end
    % Save backtracked image
    imwrite(uint16(Backtracked),[Path filesep timing '_' num2str(embryo) filesep 'FoxA2_ventral(binarized-backtracked4h).tif']);

    % Open binarized FoxA2_dorsal
    Image=imread([Path filesep timing '_' num2str(embryo) filesep 'FoxA2_dorsal.tif']);
    % Filling up the backtracked image
    Backtracked=nan(size(Image));
    for x=1:size(Image,2)
        for y=1:size(Image,1)
            X=LagrangePosX(y,x);
            Y=LagrangePosY(y,x);
            if X>=1 && X<=size(Image,2) && Y>=1 && Y<=size(Image,1) && ~isnan(X) && ~isnan(Y)
                Backtracked(y,x)=Image(Y,X);
            end
        end
    end
    % Save backtracked image
    imwrite(uint16(Backtracked),[Path filesep timing '_' num2str(embryo) filesep 'FoxA2_dorsal(backtracked4h).tif']);
    toc
end
