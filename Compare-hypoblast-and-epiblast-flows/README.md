# Compare hypoblast and epiblast flows

**Context and aims:** Chimeras made with hypoblast cells expressing tdTomato-Myosin grafted on a host expressing memGFP are imaged from the epiblast side. MemGFP signal can be used to extract epiblast flows by PIV, while tdTomato can be used to extract hypoblast flows by manual tracking. This code aims at generating average maps of epiblast and hypoblast flows to compare them (found in Figure 2B-D, Movie 4 and Supplementary Figure 3A, 3B and 3D).


## Step1: Manual tracking of hypoblast cells (tdTomato-Myosin signal) 
-	**Tool**:  Fiji (`Plugins→Tracking→Manual Tracking`)
-	**Input data**: Timelapse movie of tdTomato-Myosin signal (like in the following example: (`Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Graft.tif`)
-	**Instructions**: In Fiji, open the movie of tdTomato-Myosin signal. Track all the detectable hypoblast islands in time using the Manual Tracking plugin `Plugins→Tracking→Manual Tracking`). Save the tracking as (`Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Tracking.csv`).
-	**Output**: `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Tracking.csv`, containing hypoblast islands tracking data.

## Step2: PIV tracking of the epiblast
-	**Tool**: PIV analysis method described in *Saadaoui & al., Science (2020)*
-	**Input data**: Timelapse movie of epiblast dynamics (like in the following example:`Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Host.tif`)
-	**Ouput**: `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/data/Host-flows.h5`, containing epiblast PIV tracking data. `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/movies/Host-VEC.tif`, movie displaying velocity vectors.

## Step3: Spatial and temporal alignment
-	**Tool**: Fiji and Excel
-	**Input data**: `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/movies/Host-VEC.h5`, movie displaying velocity vectors.
-	**Instructions**: On the movie displaying velocity vectors (`Compare-hypoblast-and-epiblast-flows/Example/Chimera1/movies/Host-VEC.tif`), click on the center of counter-rotating flows at 8h. Save the coordinates as `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Registration.csv`.
In Excel, create a column containing the timing in h of each frame of the movie. Save it as `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Timing.csv`.
-	**Output**: `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Registration.csv` and `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Timing.csv`, containing respectively information for spatial and temporal registration.

## Step4: Pooling data for all animals and generate average maps
-	**Tool**: Matlab (`Compare-hypoblast-and-epiblast-flows/Code/CompareHypoblastEpiblastFlows.m`)
-	**Input**: For each movie, registration information (`Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Registration.csv` and `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Timing.csv`), hypoblast tracking data (`Compare-hypoblast-and-epiblast-flows/Example/Chimera1/Tracking.csv`) and epiblast PIV data (`Compare-hypoblast-and-epiblast-flows/Example/Chimera1/data/Host-flows.h5`)
-	**Instructions**: In Matlab, open `Compare-hypoblast-and-epiblast-flows/Code/CompareHypoblastEpiblastFlows.m`. In the “Parameters” section, set *Path* so that it corresponds to the path finishing with `PIV-analysis_Compare-with-hypoblast-flows/Example`, and *PathData* so that it corresponds to the path finishing with `Compare-hypoblast-and-epiblast-flows/Data`. Run the code.
-	**Description of the code**: `CompareHypoblastEpiblastFlows.m` first loads epiblast and hypoblast tracking data, as well as spatial and temporal alignment information for each movie. For each hypoblast island, located at (x,y) at 4h, it will store the displacement from 4h to 12h at the location (x,y) in the epiblast and the hypoblast. It will as well compute at each timepoint from 4h to 12h the differential motion between the hypoblast and the epiblast of the given hypoblast island and sum it up to get differential motion from 4h to 12h.
`CompareHypoblastEpiblastFlows.m` will then register in space the data of all the movies, to pool them and average them among all the animals. It then generates average velocity field from 4h to 12h in the epiblast and the hypoblast, as well as an average map of differential motion (Figure 2C-C’). A grid representation of the velocity fields is also plotted (Figure 2D-D’).
For each chimera, `CompareHypoblastEpiblastFlows.m` also generates .txt files required to plot the tracks of the epiblast, hypoblast, and differential motion (used in Supplementary Figure 3A, 3B and 3D, see Step5).
-	**Output**: Windows1-3: all velocity vectors of epiblast, hypoblast or differential motion from all the registered movies, pooled together. Windows 4-6: average velocity fields of epiblast, hypoblast and differential motion. Windows 7-9: associated deformation grids.
`Compare-hypoblast-and-epiblast-flows/Example/Chimera1/values`: folder containing .txt files used for plotting epiblast, hypoblast and differential motion tracks (used in **Step5**).

## Step5: For a given movie, plot epiblast, hypoblast and differential motion track (to generate an illustration of the dataset, like in Movie 4 or Supplementary Figure 3A, 3B and 3D)
-	**Tool**: Fiji (`Compare-hypoblast-and-epiblast-flows/Code/PlotTrack.ijm`)
-	**Input data**: `Compare-hypoblast-and-epiblast-flows/Example/Chimera1/values` content (.txt files)
-	**Instructions**: Run `Compare-hypoblast-and-epiblast-flows/Code/PlotTrack.ijm` in Fiji. Indicate PIV-analysis_Compare-with-hypoblast-flows/Example/Chimera1/values when a path is asked. The code will generate the illustrations of the different tracks (epiblast, hypoblast and differential tracks and dots).
-	**Output**: Illustration movies of epiblast, hypoblast and differential motion tracks for a given movie


>Requirements: 
>Use Matlab (2023 version), and Excel (or in any similar software able to generate .csv files).
