%% ArchetypePattern
% Collect and average cropped circular patterns

clearvars
close all

%% Parameters
% Path where to find embryo subfolders
Path='\\gaia.pasteur.fr\MVS_DATA1\Aurelien\PAPER\23-02-24_Hypoblast-paper\CodeAvailability\Staining-signal-quantification\Data';
% Timings to consider
Timings={'2h' '4h' '6h' '8h'};
% Indexes of the embryos associated with each timing
Indexes={1:7 1:6 1:11 1:7};
% Size of the images
Size=600;

%% Binarized NODAL patterns
close all
for timing=1:4
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
    colormap jet
    colorbar
    title([Timings{timing} ' ventral'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[Path filesep 'NODAL-ventral' num2str(Timings{timing}) '.png'],'-r300');

    hfig=figure();
    imagesc(mean(PoolNODALdorsal,3))
    colormap jet
    colorbar
    title([Timings{timing} ' dorsal'])
    caxis([0 1])
    daspect([1 1 1])
    xticks([]);
    yticks([]);
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    print(hfig,'-dpng',[Path filesep 'NODAL-dorsal' num2str(Timings{timing}) '.png'],'-r300');

end
