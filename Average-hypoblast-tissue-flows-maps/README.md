# **Average hypoblast tissue flows maps**

**Aim: Generate average maps of hypoblast flows by averaging PIV data from many animals.** 


## Step1: Generate the metadata for space and time alignment for each movie 
-	**Tool**:  Fiji and Excel
-	**Input data**: Timelapse movie of hypoblast dynamics (like in the following example: `Average-hypoblast-tissue-flows-maps/Example/Embryo1/MAX.tif`), and associated visualization of vector fields (`Average-hypoblast-tissue-flows-maps/Example/Embryo1/movies/MAX-VEC.tif`)
-	**Instructions**: In Fiji, open the movie displaying vector fields (`Average-hypoblast-tissue-flows-maps/Example/Embryo1/movies/MAX-VEC.tif`). Display the frame corresponding to 8 h and measure the coordinates of the two centers of the counter-rotating flows. Save them in `Average-hypoblast-tissue-flows-maps/Example/Embryo1/CoordinatesCentersRotation.csv` (as in the example).
In Excel, generate a column corresponding to the timing post-laying (in hour) of each frame of the movie. Save it in `Average-hypoblast-tissue-flows-maps/Example/Embryo1/Timing.csv`.
-	**Output**: `CoordinatesCentersRotation.csv` and `Timing.csv`, containing respectively information relative to spatial and temporal alignment.

## Step2: Generate a mask identifying debris masking the signal and preventing correct PIV calculation
-	**Tool**: Fiji
-	**Input data**: Timelapse movie of hypoblast dynamics (`Average-hypoblast-tissue-flows-maps/Example/Embryo1/MAX.tif`)
-	**Instructions**: Open the movie in Fiji. Generate a mask identifying for each time point the position of the debris that prevent PIV calculation. Debris can be identified by their low fluorescence signal by applying a threshold (Image→Adjust→Threshold) after applying a median filter to the image (Process→Filters→Median) or by manually defining the contour of the debris (Freehand Selection, then Edit→Fill). Save the mask as `Average-hypoblast-tissue-flows-maps/Example/Embryo1/Mask.tif`.
-	**Ouput**: `Average-hypoblast-tissue-flows-maps/Example/Embryo1/Mask.tif`, a binarized version of the movie where regions to filter out (masked by debris) are 1.

## Step3: Pool PIV data and generate average maps
-	**Tool**: Matlab (`Average-hypoblast-tissue-flows-maps/Code/HypoblastFlowAnalyzer.m`)
-	**Input data**: Timelapse movie of hypoblast dynamics (`Average-hypoblast-tissue-flows-maps/Example/Embryo1/MAX.tif`), associated PIV analysis tracking data using PIV analysis pipeline developed in *Saadaoui & al., Science (2020)* (`Average-hypoblast-tissue-flows-maps/Example/Embryo1/data/MAX-flows.h5`), mask identifying debris (`Average-hypoblast-tissue-flows-maps/Example/Embryo1/Mask.tif`), spatio-temporal alignment data (`Average-hypoblast-tissue-flows-maps/Example/Embryo1/CoordinatesCentersRotation.csv` and `Average-hypoblast-tissue-flows-maps/Example/Embryo1/Timing.csv`)
-	**Instructions**: In Matlab, open `Average-hypoblast-tissue-flows-maps/Code/HypoblastFlowAnalyzer.m`. In the “Parameters” section, change *Path* to make it correspond to the path ending with `Average-hypoblast-tissue-flows-maps/Example`, and *PathOut* to make it correspond to the path ending with `Average-hypoblast-tissue-flows-maps/Data`. Run the code.
-	**Description of the code**: `HypoblastFlowAnalyzer.m` first uploads PIV data for all the movies (here, only one example is provided), interpolates the data in (x,y,t) to align all the movies together and save the pooled dataset. `HypoblastFlowAnalyzer.m` then reads the pooled dataset (here provided for all the animals) to generate average maps (like in Figure 1 and Movie 1). Only (x,y) regions covered in at least 3 movies are displayed.
-	**Preset parameters**: *tmin*, *tmax*, *tstep*, *Xmin*, *Xmax*, *Xstep*, *Ymin*, *Ymax*, *Ystep*: parameters related to the size and spacing of the interpolation grid (x,y,t)
-	**Ouput**: Average maps (Figure 1, Movie 1)


>Requirements: 
>Use Matlab (2023 version), and Excel (or in any similar software able to generate .csv files).
