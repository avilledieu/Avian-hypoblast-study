# Staining signal backtracking

**Context and aims:** For each imaged and stained embryo, the first aim is to backtrack 8h-NODAL RNAs level pattern from 8h to 4h, using flows from 4h to 8h quantified by PIV. Then, the final aim is to generate an archetype map of 8h-NODAL RNAs levels backtracked to 4h by averaging the signals from all the embryos (Supplementary Fig. 6a).

**Input data:** 
- For each embryo:
  -	a timelapse movie of epiblast morphogenesis stopping at 8h, analyzed by PIV
  -	associated epiblast and hypoblast images of NODAL RNAs staining at 8h (see `Staining-signal-quantification` section).
- `DiffEpiHypo.mat` matrix, which contains the differential motion between hypoblast and epiblast from 4h to 8h, calculated in the `PIV-analysis/Compare-with-hypoblast-flows` section.



## Step1: For each embryo, align the staining images with the last frame of the timelapse movie 
-	**Tool**:  Fiji
-	**Input data**: Last frame (corresponding to 8h) of the timelapse movie (`Staining-signal-backtracking/Example/8h_1/Live/MAX.tif`), images of the staining (`Staining-signal-backtracking/Example/8h_1/NODAL_dorsal.tif` and `NODAL_ventral.tif`).
-	**Output**: Staining images aligned with the last frame of the timelapse movie.
-	**Instructions**: **(1) Flip horizontally the hypoblast view of the staining**: Hypoblast view is reversed as compared with epiblast timelapse movie view. Use `Image→Transform→Flip Horizontally` to correct for the orientation of the hypoblast staining image. **(2) Rescaling**: Using 2 recognizable landmarks, measure a representative length of the embryo in the last frame of the timelapse and in the staining (make sure `Pixel size` is set to 1 in `Image→Properties`). Calculate the ratio between the two values and use this ratio as a parameter in `Image→Scale`, to rescale the staining image to the size of the timelapse movie. **(3) Rotation and translation**: using recognizable landmarks, calculate the tilt of the staining image as regards the tilt of the last frame of the timelapse movie, and correct for it using `Image→Transform→Rotate`. Similarly, calculate the translation to apply to align the staining image with the last frame of the timelapse movie using `Image→Transform→Translate`. **(4) Saving**: Save the aligned staining image by overwriting the previous unaligned images.

## Step2: For each embryo, binarize the staining images
-	**Tool**: Fiji
-	**Input data**: Aligned images of the staining (`Staining-signal-backtracking/Example/8h_1/NODAL_dorsal.tif` and `NODAL_ventral.tif`).
-	**Ouput**: `Staining-signal-backtracking /Example/8h_1/NODAL_dorsal(binarized).tif` and `NODAL_ventral(binarized).tif`, which are binarized version of NODAL signal.
-	**Instructions**: Open staining aligned image. Zoom in a region where isolated dots corresponding to NODAL positive signal can be seen. Apply a threshold, manually adjusted (`Image→Adjust→Threshold`) so that isolated dots are well segmented (such thresholding may result in a saturation in regions containing a lot of overlapped dots, especially for stages when NODAL expression is high). Save the binarized image (this procedure has already been done for generating archetypal NODAL RNAs levels in `Staining-signal-quantification section`).

## Step3: For each embryo, spatial-temporal alignment of the movie (identical to Step1 in `Archetypal-PIV-maps` section)
-	**Tool**: Fiji
-	**Input data**: Timelapse movie of hypoblast dynamics (like in the following example: `Staining-signal-backtracking/Example/8h_1/Live/MAX.tif`), and associated visualization of vector fields (`Staining-signal-backtracking/Example/8h_1/Live/movies/MAX-VEC.tif`)
-	**Ouput**: `CoordinatesCentersRotation.csv` and `Timing.csv`, containing respectively information relative to spatial and temporal alignment.
-	**Instructions**: In Fiji, open the movie displaying vector fields (`Staining-signal-backtracking/Example/8h_1/Live/MAX-VEC.tif`). Display the last frame corresponding to 8h and measure the coordinates of the two centers of the counter-rotating flows. Save them in `Staining-signal-backtracking/Example/8h_1/Live/CoordinatesCentersRotation.csv` (as in the example).

## Step4: For each embryo, backtraking of the 8h-NODAL RNAs pattern (epiblast and hypoblast ones) to 4h using PIV calculated flows
-	**Tool**: Matlab (`Staining-signal-backtracking/Code/BacktringPattern.m`)
-	**Input data**: PIV calculated flows (`Staining-signal-backtracking/Example/8h_1/Live/data/MAX-flows.h5`), binarized aligned NODAL patterns in the epiblast and the hypoblast at 8h (`Staining-signal-backtracking/Example/8h_1/NODAL_dorsal(binarized).tif` and `NODAL_ventral(binarized).tif`), spatio-temporal alignment metadata (`CoordinatesCentersRotation.csv` and `Timing.csv`), and `DiffEpiHypo.mat` matrix (which contains the differential motion between hypoblast and epiblast from 4h to 8h, calculated in the `PIV-analysis_Compare-with-hypoblast-flows` section).
-	**Ouput**: For each embryo, 8h-NODAL RNAs levels backtracked to 4h in the epiblast and hypoblast sides (`Staining-signal-backtracking/Example/8h_1/NODAL_dorsal(binarized-backtracked4h).tif` and `NODAL_ventral(binarized-backtracked4h).tif`)
-	**Instructions**: Open `BacktrackingPattern.m` in Matlab. Specify the `Path` where to find the subfolders named ‘8h_1’, ‘8h_2’,… containing the data for each embryo and run the program. For each embryo, the program will calculate epiblast motion from 4h to 8h. It will add the average differential motion between hypoblast and epiblast contained in `DiffEpiHypo.mat` to epiblast motion (aligning it using `CoordinatesCentersRotation.csv`) for estimating the motion of the hypoblast. It well then backtrack 8h patterns to 4h and save the backtracked epiblast and hypoblast NODAL levels patterns for each embryo.
  
## Step5: For each embryo, spatial alignment of the backtracked staining (identical to Step 1 of `Staining-signal-quantification` section)
-	**Tool**: Fiji
-	**Input data**: 2h frame of the timelapse movie aligned with the staining (in `Staining-signal-backtracking/Example/8h_1/Live/MAX.tif`)
-	**Ouput**: `Staining-signal-backtracking/Example/8h_1/Spatial-alignment.zip`, containing the 2h contours of the blastoderm and the circle extracted from it used to align staining data. `Staining-signal-backtracking/Example/8h_1/Angle.csv`, containing the tilt of the staining.
-	**Instructions**: **(1) Segment blastoderm contour at 2h**: In Fiji, open the 2h picture of the embryo, and extract the contour of the blastoderm by applying an average filter (`Process→Filters→Mean`) and thresholding the image (`Image→Adjust→Threshold`). Save the segmented blastoderm image in `Staining-signal-backtracking/Example/8h_1/MaskEmbryo.tif` and add the segmented blastoderm contour to the ROI Manager (`Wand tool`, `Edit→Selection→Add to Manager`). **(2) Extract the alignment circle from blastoderm contour**: Fit a circle to the blastoderm contour (`Edit→Selection→Fit Circle`). Measure its radius and decrease it by 10% (`Edit→Selection→Enlarge`) to obtain the circle used for spatial alignment. Add it to the ROI Manager and save it as `Staining-signal-backtracking/Example/8h_1/Spatial-alignment.zip`.
**(3) Measure the tilt of the staining**: In Fiji, measure the tilt of the staining, so that NODAL-positive crescent is left-right symmetrical. Save the measured tilt angle in `Staining-signal-backtracking/Example/8h_1/Angle.csv`.

## Step6: For each embryo, crop, resize and rotate the backtracked patterns (identical to Step 3 of `Staining-signal-quantification` section)
-	**Tool**: Fiji
-	**Input data**: Binarized backtracked NODAL levels images (`Staining-signal-backtracking/Example/8h_1/NODAL_dorsal(binarized).tif` and `NODAL_ventral(binarized).tif`). Spatial alignment information (`Staining-signal-backtracking/Example/8h_1/Spatial-alignment.zip` and `Angle.csv`)
-	**Ouput**: Aligned and resized staining images `Staining-signal-backtracking/Example/8h_1/NODAL_dorsal(binarized-crop-rotated).tif` and `NODAL_ventral(binarized-crop-rotated).tif`
-	**Instructions**: Open the binarized images in Fiji. Open the ROIs contained in `Spatial-alignment.zip`, and crop the images using the alignment circle (`Image→Crop`). Resize the image so that it is 600X600 pixels (`Image→Scale`). Rotate the cropped and donwscaled image using the measured angle in `Angle.csv` (`Image→Transform→Rotate`). Save the cropped and rotated image.

## Step7: Pool the backtracked cropped binarized patterns of all embryos and average them (identical to Step 4 of `Staining-signal-quantification` section)
-	**Tool**: Matlab (`Staining-signal-backtracking/Code/ArchetypePattern.m`)
-	**Input data**: `Staining-signal-backtracking/Example/8h_1/NODAL_dorsal(binarized-crop-rotated).tif` and `NODAL_ventral(binarized-crop-rotated).tif`, for several embryos and several timings
-	**Ouput**: Average map of 8h-NODAL mRNA localization in the epiblast and the hypoblast backtracked to 4h (Supplementary Fig. 6a).
-	**Instructions**: In Matlab, open `Staining-signal-backtracking /Code/ArchetypePattern.m`. Adjust `Path` so that it corresponds to the path ending with `Staining-signal-backtracking/Data`. Run the code.

>Requirements: 
>Download Fiji (https://imagej.net/software/fiji/downloads)
>
>Use Matlab (2023 version)
