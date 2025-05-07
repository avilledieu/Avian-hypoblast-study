# Staining signal quantification

**Aim:** Generate an archetypal map of *NODAL* mRNA localization, by averaging the signals from many different animals stained at a given timing (Figure 3B).


## Step1: Spatial alignment of the staining
-	**Tool**:  Fiji
-	**Input data**: Hoechst signal (nuclear signal) of the staining of the embryo at 2h (`Staining-signal-quantification/Example/2h_1/HOECHST_ventral.tif`), or 2h frame of the timelapse movie aligned with the staining (memGFP signal).
-	**Instructions**: (1) Segment blastoderm contour at 2h: In Fiji, open the 2h picture of the embryo, and extract the contour of the blastoderm by applying an average filter (`Process→Filters→Mean`) and thresholding the image (`Image→Adjust→Threshold`). Save the segmented blastoderm image in `Staining-signal-quantification/Example/2h_1/MaskEmbryo.tif` and add the segmented blastoderm contour to the ROI Manager (`Wand tool`, `Edit→Selection→Add to Manager`).
(2) Extract the alignment circle from blastoderm contour: Fit a circle to the blastoderm contour (`Edit→Selection→Fit Circle`). Measure its radius and decrease it by 10% (Edit→Selection→Enlarge) to obtain the circle used for spatial alignment. Add it to the ROI Manager and save it as `Staining-signal-quantification/Example/2h_1/Spatial-alignment.zip`.
(3) Measure the tilt of the staining: In Fiji, measure the tilt of the staining, so that NODAL-positive crescent is left-right symmetrical. Save the measured tilt angle in `Staining-signal-quantification/Example/2h_1/Angle.csv`.
-	**Output**: `Staining-signal-quantification/Example/2h_1/Spatial-alignment.zip`, containing the 2h contours of the blastoderm and the circle extracted from it used to align staining data. `Staining-signal-quantification/Example/2h_1/Angle.csv`, containing the tilt of the staining.

## Step2: Binarization of NODAL HCR-RNA-FISH signal
-	**Tool**: Fiji
-	**Input data**: *NODAL* HCR-RNA-FISH projected signal for the epiblast (`Staining-signal-quantification/Example/2h_1/NODAL_dorsal.tif`) and the hypoblast side (`Staining-signal-quantification/Example/2h_1/NODAL_ventral.tif`).
-	**Instructions**:  Open *NODAL* projected signal. Zoom in a region where isolated dots corresponding to *NODAL* positive signal can be seen. Apply a threshold, manually adjusted (`Image→Adjust→Threshold`) so that isolated dots are well segmented (such thresholding may result in a saturation in regions containing a lot of overlapped dots, especially for stages when NODAL expression is high). Save the binarized image.
-	**Output**: `Staining-signal-quantification/Example/2h_1/NODAL_dorsal(binarized).tif` and `Staining-signal-quantification/Example/2h_1/NODAL_ventral(binarized).tif`, which are binarized version of NODAL signal, which can subsequently be used for averaging and comparing different timings.

## Step3: Crop, resize and rotate staining using alignment circle
-	**Tool**: Fiji
-	**Input data**: Binarized *NODAL* HCR images (`Staining-signal-quantification/Example/2h_1/NODAL_dorsal(binarized).tif` and `NODAL_ventral(binarized).tif`).
Spatial alignment information (`Staining-signal-quantification/Example/2h_1/Spatial-alignment.zip` and `Angle.csv`)
-	**Instructions**: Open the binarized images in Fiji. Open the ROIs contained in `Spatial-alignment.zip`, and crop the images using the alignment circle (`Image→Crop`). Resize the image so that it is 600X600 pixels (Image→Scale). Rotate the cropped and donwscaled image using the measured angle in Angle.csv (`Image→Transform→Rotate`). Save the cropped and rotated image.
-	**Output**: `Staining-signal-quantification/Example/2h_1/NODAL_dorsal(binarized-crop-rotated).tif` and `NODAL_ventral(binarized-crop-rotated).tif`

## Step 4: Pool the cropped binarized images of the staining and average them
-	**Tool**: Matlab (`Staining-signal-quantification/Code/ArchetypePattern.m`)
-	**Input data**: `Staining-signal-quantification/Data/2h_1/ NODAL_dorsal(binarized-crop-rotated).tif` and `NODAL_ventral(binarized-crop-rotated).tif`, for several embryos and several timings
-	**Instructions**: In Matlab, open `Staining-signal-quantification/Code/ArchetypePattern.m`. Adjust *Path* so that it corresponds to the path ending with `Staining-signal-quantification/Data`. Run the code.
-	**Output data**: Average map of *NODAL* mRNA localization in the epiblast and the hypoblast at 2, 4, 6 and 8h (Figure 3B)

>Requirements: 
>Download Fiji (https://imagej.net/software/fiji/downloads)
>Use Matlab (2023 version)
