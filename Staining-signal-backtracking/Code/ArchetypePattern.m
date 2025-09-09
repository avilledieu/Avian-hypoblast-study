%% ArchetypePattern
% Collect and average cropped circular patterns (early embryo, early
% pattern or pattern backtracked to early embryo).

clearvars
close all

%% Parameters
% Path where to find embryo subfolders
Path='\\gaia.pasteur.fr\MVS_DATA1\Aurelien\9_Fixed-samples\24-06-06_SUMMARY-Nodal-FoxA2-patterns-in-time';
% Timings to consider
Timings={'8h'};
% Indexes of the embryos associated with each timing
Indexes={1:7};
% Indexes of the embryos associated with each timing (FoxA2)
IndexesFoxA2={1:6};
% Size of the images
Size=600;
% Average filter radius used for FoxA2 averaging (in pixel on the resized image)
AverageRadius=20;
% Output folder
PathOut='\\gaia.pasteur.fr\MVS_DATA1\Aurelien\9_Fixed-samples\24-06-06_SUMMARY-Nodal-FoxA2-patterns-in-time\PLOTS';

%% Binarized NODAL patterns
close all
for timing=1
    Index=Indexes{timing};
    %% Collect all images
    % Initialization of the pooling matrices
    PoolNODALventral=nan(Size,Size,length(Index));
    PoolNODALdorsal=nan(Size,Size,length(Index));
    for embryo=Index
        % NODAL_ventral
        temp=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'NODAL_ventral(binarized-crop-rotated).tif']);
        temp=double(temp);
        %temp=imboxfilt(temp,FilterSize);
        PoolNODALventral(:,:,embryo)=temp./255;

        % NODAL_dorsal
        temp=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'NODAL_dorsal(binarized-crop-rotated).tif']);
        temp=double(temp);
        %temp=imboxfilt(temp,FilterSize);
        PoolNODALdorsal(:,:,embryo)=temp./255;
    end

    hfig=figure();
    imagesc(mean(PoolNODALventral,3))
    colormap gray
    colorbar
    title([Timings{timing} ' ventral'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'NODAL-ventral' num2str(Timings{timing}) '.png'],'-r300');

    hfig=figure();
    imagesc(mean(PoolNODALdorsal,3))
    colormap gray
    colorbar
    title([Timings{timing} ' dorsal'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'NODAL-dorsal' num2str(Timings{timing}) '.png'],'-r300');

    % for animal=Index
    % figure(animal)
    % imagesc(PoolNODALdorsal(:,:,animal))
    % colormap jet
    % colorbar
    % caxis([0 6])
    % end
end


%% FoxA2 patterns
close all
for timing=1:4
    Index=IndexesFoxA2{timing};
    %% Collect all images
    % Initialization of the pooling matrices
    PoolFoxA2ventral=nan(Size,Size,length(Index));
    PoolFoxA2dorsal=nan(Size,Size,length(Index));
    for embryo=Index
        % FoxA2_ventral averaged signal
        Image=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'FoxA2_ventral(crop-averaged-rotated).tif']);
        Image=double(Image);
        % FoxA2_ventral binary map, to reintroduce NaN values
        Mask=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'FoxA2_ventral(binarized-crop-rotated).tif']);
        Mask=imboxfilt(Mask,AverageRadius+1);
        Mask=Mask(Mask>20);
        % Apply Mask
        Image(Mask==0)=nan;
        % Normalize each image so that the maximal value corresponds to 1
        Image=Image./max(Image(:));
        % Add to the pooling matrix
        PoolFoxA2ventral(:,:,embryo)=Image;


        % FoxA2_dorsal
        Image=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'FoxA2_dorsal(crop-rotated).tif']);
        Image=double(Image);
        % Smooth the image
        Image=imboxfilt(Image,AverageRadius+1);
        % Load noise and subtract it
        Noise=csvread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'Noise_FoxA2.csv']);
        Image=Image-Noise;
        Image(Image<0)=0;
        % Add to pooling matrix
        PoolFoxA2dorsal(:,:,embryo)=Image;
    end

    hfig=figure();
    imagesc(mean(PoolFoxA2ventral,3,'omitnan'))
    colormap gray
    colorbar
    title([Timings{timing} ' ventral'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'FoxA2-ventral' num2str(Timings{timing}) '.png'],'-r300');

    hfig=figure();
    imagesc(mean(PoolFoxA2dorsal,3))
    colormap gray
    colorbar
    title([Timings{timing} ' dorsal'])
    daspect([1 1 1])
    caxis([0 15])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'FoxA2-dorsal' num2str(Timings{timing}) '.png'],'-r300');

    % for animal=Index
    % figure(animal)
    % imagesc(PoolNODALdorsal(:,:,animal))
    % colormap jet
    % colorbar
    % caxis([0 6])
    % end
end


%% Backtracked FoxA2 patterns
close all
for timing=1
    Index=IndexesFoxA2{timing};
    %% Collect all images
    % Initialization of the pooling matrices
    PoolFoxA2ventral=nan(Size,Size,length(Index));
    PoolFoxA2dorsal=nan(Size,Size,length(Index));
    for embryo=Index
        % FoxA2_ventral averaged signal
        Image=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'FoxA2_ventral(backtracked-dorsal-crop-rotated-averaged).tif']);
        Image=double(Image);
        % FoxA2_ventral binary map, to reintroduce NaN values
        Mask=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'FoxA2_ventral(binarized-backtracked-crop-rotated).tif']);
        Mask=imboxfilt(Mask,AverageRadius+1);
        Mask=Mask(Mask>20);
        % Apply Mask
        Image(Mask==0)=nan;
        % Normalize each image so that the maximal value corresponds to 1
        Image=Image./max(Image(:));
        % Add to the pooling matrix
        PoolFoxA2ventral(:,:,embryo)=Image;


        % FoxA2_dorsal
        Image=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'FoxA2_dorsal(backtracked-crop-rotated).tif']);
        Image=double(Image);
        % Smooth the image
        Image=imboxfilt(Image,AverageRadius+1);
        % Load noise and subtract it
        Noise=csvread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'Noise_FoxA2.csv']);
        Image=Image-Noise;
        Image(Image<0)=0;
        % Add to pooling matrix
        PoolFoxA2dorsal(:,:,embryo)=Image;
    end

    hfig=figure();
    imagesc(mean(PoolFoxA2ventral,3,'omitnan'))
    colormap gray
    colorbar
    title([Timings{timing} ' ventral'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'FoxA2-ventral' num2str(Timings{timing}) '(backtracked).png'],'-r300');

    hfig=figure();
    imagesc(mean(PoolFoxA2dorsal,3))
    colormap gray
    colorbar
    title([Timings{timing} ' dorsal'])
    daspect([1 1 1])
    caxis([0 15])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'FoxA2-dorsal' num2str(Timings{timing}) '(backtracked).png'],'-r300');
end


%% Binarized backtracked dorsal NODAL patterns
close all
for timing=1
    Index=Indexes{timing};
    %% Collect all images
    % Initialization of the pooling matrices
    PoolNODALventral=nan(Size,Size,length(Index));
    PoolNODALdorsal=nan(Size,Size,length(Index));
    for embryo=Index
        % NODAL_dorsal
        temp=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'NODAL_dorsal(binarized-backtracked4h-crop-rotated).tif']);
        temp=double(temp);
        %temp=imboxfilt(temp,FilterSize);
        PoolNODALdorsal(:,:,embryo)=temp./255;
    end

    hfig=figure();
    imagesc(mean(PoolNODALdorsal,3))
    colormap gray
    colorbar
    title([Timings{timing} ' dorsal'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'NODAL-dorsal' num2str(Timings{timing}) '(backtracked4h).png'],'-r300');

    for animal=1:size(PoolNODALdorsal,3)
        figure(animal)
        imagesc(PoolNODALdorsal(:,:,animal))
        colormap jet
        colorbar
        title([Timings{timing} ' ventral'])
        caxis([0 1])
        daspect([1 1 1])
        xticks([]);
        yticks([]);
        set(findall(gcf,'-property','FontSize'),'FontSize',13)
    end

    hfig=figure();
    imagesc(mean(PoolNODALdorsal(:,:,[1 3:7]),3))
    colormap gray
    colorbar
    title([Timings{timing} ' dorsal'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'NODAL-dorsal' num2str(Timings{timing}) '(backtracked)-without2.png'],'-r300');
end


%% Binarized backtracked ventral NODAL patterns
close all
for timing=1
    Index=Indexes{timing};
    %% Collect all images
    % Initialization of the pooling matrices
    PoolNODALventral=nan(Size,Size,length(Index));
    PoolNODALdorsal=nan(Size,Size,length(Index));
    for embryo=Index
        % NODAL_ventral
        temp=imread([Path filesep Timings{timing} '_' num2str(embryo) filesep 'NODAL_ventral(binarized-backtracked4h-crop-rotated).tif']);
        temp=double(temp);
        %temp=imboxfilt(temp,FilterSize);
        PoolNODALventral(:,:,embryo)=temp./255;
    end

    hfig=figure();
    imagesc(mean(PoolNODALventral,3))
    colormap gray
    colorbar
    title([Timings{timing} ' ventral'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'NODAL-ventral' num2str(Timings{timing}) '(backtracked).png'],'-r300');

    for animal=1:size(PoolNODALventral,3)
        figure(animal)
        imagesc(PoolNODALventral(:,:,animal))
        colormap jet
        colorbar
        title([Timings{timing} ' ventral'])
        caxis([0 1])
        daspect([1 1 1])
        xticks([]);
        yticks([]);
        set(findall(gcf,'-property','FontSize'),'FontSize',13)
    end

    hfig=figure();
    imagesc(mean(PoolNODALventral(:,:,[1 3:7]),3))
    colormap gray
    colorbar
    title([Timings{timing} ' ventral'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[PathOut filesep 'NODAL-ventral' num2str(Timings{timing}) '(backtracked)-without2.png'],'-r300');

end

