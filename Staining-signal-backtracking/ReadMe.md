# Staining signal backtracking

**Context and aims:** For each imaged and stained embryo, backtrack 8h-NODAL RNAs level pattern from 8h to 4h, using flows from 4h to 8h quantified by PIV. Generate an archetype map of 8h-NODAL RNAs levels backtracked to 4h by averaging the signals from all the embryos (Supplementary Fig. 6a).

**Input data:** 
- For each embryo:
  -	a timelapse movie of epiblast morphogenesis stopping at 8h, analyzed by PIV
  -	associated epiblast and hypoblast images of NODAL RNAs staining at 8h (see Staining-signal-quantification section).
- `DiffEpiHypo.mat` matrix, which contains the differential motion between hypoblast and epiblast from 4h to 8h, calculated in the `PIV-analysis/Compare-with-hypoblast-flows` section.



## Step1: For each embryo, align the staining images with the last frame of the timelapse movie 
-	**Tool**:  Fiji
-	**Input data**: Last frame (corresponding to 8h) of the timelapse movie (Staining-signal-backtracking/Example/8h_1/Live/MAX.tif), images of the staining (Staining-signal-backtracking/Example/8h_1/NODAL_dorsal.tif and NODAL_ventral.tif).
-	**Output**: Downscaled and binarized image of the Hoechst channel (nuclear marker), saved in “TopoMapDetection” folder. This image is going to be used in the following step of the pipeline to detect the most superficial plane in the hyperstack.
-	**Instructions**: In Fiji, run the code `Staining-signal-projection/Code/1_Downscale-and-binarize.ijm`. When asked in Fiji, provide the path leading to the `Staining-signal-projection/Example`. Files ending with “_HOECHST.tif” contained in this folder will be open and downscaled. The user is asked to threshold the downscaled image, so that the most superficial z-plane positive for Hoechst for each pixel (x,y) is the first white pixel when travelling from the top to the bottom of the hyperstack.
-	**Preset parameters**: *RollingBallRadius* and *MeanFilter* (parameters associated to preprocessing filters to reduce noise), *ScalingFactor* (downscaling factor).

## Step2: Projection of the signals using the binary degraded image
-	**Tool**: Matlab (`Staining-signal-projection/Code/ProjectionSuperficialSurface.m`)
-	**Input data**: Hyperstacks to project (contained in `Staining-signal-extraction/Example`) and associated downgraded and binarized images (contained in `Staining-signal-projection/Example/TopoMapDetection`).
-	**Ouput**: Projection of each channel around the most supercial z-plane, contained in `Staining-signal-projection/Example/Projected`.
-	**Instructions**:  In Matlab, open `Staining-signal-projection/Code/ProjectionSuperficialSurface.m`. In the “Parameters” section, change *Path* so that it corresponds to the path ending with `Staining-signal-projection/Example`. Run the code.
-	**Parameters**: *ApicalShift* and *BasalShift* (number of z-plane around which the signal is projected, preset to +1 around the detected z-plane)

## Step3: Alignment of epiblast and hypoblast signals
-	**Tool**: Fiji
-	**Input data**: Epiblast and hypoblast side projections (found in `Staining-signal-projection/Example/Projected`)
-	**Ouput**: Aligned version of epiblast and hypoblast side projections. An example can be found in `Staining-signal-projection/Example/Projection/Aligned`.
-	**Instructions**: Horizontally flip hypoblast projection (`Image→Transform→Flip Horizontally`). By identifying landmarks on epiblast and hypoblast images, rotate (`Image→Transform→Rotate`) and translate (`Image→Transform→Translate`) hypoblast image to align it to epiblast image. 
Using a similar method, projections of staining signals can be aligned to the last frame of the live timelapse movie of the same embryo used for staining (used in Figure 3B or in Supplementary Figure 5).

>Requirements: 
>Download Fiji (https://imagej.net/software/fiji/downloads)
>Use Matlab (2023 version)
