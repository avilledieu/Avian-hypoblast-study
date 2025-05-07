// PlotTrack
// Plot the tracks gathered and calculated using DifferenceSpeedManualTracking.m

// Parameters ///////////////////////////////////////////////////////////////////////////////
// Diameter of the tracked points (in pixels)
Diameter=5;
/////////////////////////////////////////////////////////////////////////////////////////////


// Code /////////////////////////////////////////////////////////////////////////////////////
// Directory selection
path=getDirectory("Select 'values' directory (contains .txt output files)");

// Output directory creation
pathOut=path+"plots"+File.separator;
File.makeDirectory(pathOut);

// Activate batch mode
setBatchMode(true);

// Read values
filestring=File.openAsString(path+"Slicen.txt");
Slicen=split(filestring,",");
filestring=File.openAsString(path+"Trackn.txt");
Trackn=split(filestring,",");
filestring=File.openAsString(path+"XsDiff.txt");
XsDiff=split(filestring,",");
filestring=File.openAsString(path+"YsDiff.txt");
YsDiff=split(filestring,",");
filestring=File.openAsString(path+"XsSmoothed.txt");
XsSmoothed=split(filestring,",");
filestring=File.openAsString(path+"YsSmoothed.txt");
YsSmoothed=split(filestring,",");
filestring=File.openAsString(path+"XsEpi.txt");
XsEpi=split(filestring,",");
filestring=File.openAsString(path+"YsEpi.txt");
YsEpi=split(filestring,",");
filestring=File.openAsString(path+"ImageSize.txt");
ImageSize=split(filestring,",");
MaxFrame=File.openAsString(path+"MaxFrame.txt");
NumberTracks=File.openAsString(path+"NumberTracks.txt");


// Plot cumulative hypoblast tracks ////////////////////////////////////////////////////////////
// Generate blank image
newImage("HyperStack", "8-bit color-mode",ImageSize[0],ImageSize[1], 1,MaxFrame, 1);
ID=getImageID();
// Draw track after track
for (track=1;track<=NumberTracks;track++){
	count=0;
	for (i=0;i<Slicen.length;i++){
		if (Trackn[i]==track){
			// For the first tracked point, don't do anything
			count=count+1;
			if (count==1){
				X=XsSmoothed[i];
				Y=YsSmoothed[i];
			}
			if (count!=1){
				// For the 2nd to the last tracked point, draw the (n-1)-n segment
				X2=XsSmoothed[i];
				Y2=YsSmoothed[i];
				// Identify the slice and get to the slice
				slice=Slicen[i];
				// 
				for (s=slice;s<=MaxFrame;s++){
					selectImage(ID);
					setSlice(s);
					// Draw the segment
					makeLine(X,Y,X2,Y2);
					run("Draw", "slice");
				}
				// Keep the edge of the segment for drawing next segment
				X=X2;
				Y=Y2;
			}
		}
	}
}
// Save cumulated tracks
selectImage(ID);
saveAs("tiff",pathOut+"TracksHypo");
close();


// Plot tracked hypoblast dots  /////////////////////////////////////////////////////////////////////
// Generate blank image
newImage("HyperStack", "8-bit color-mode",ImageSize[0],ImageSize[1], 1,MaxFrame, 1);
ID=getImageID();
// Draw tracked point after tracked point
for (point=0;point<Trackn.length;point++){
	selectImage(ID);
	slice=Slicen[point];
	setSlice(slice);
	x=parseFloat(XsSmoothed[point])-Diameter/2;
	y=parseFloat(YsSmoothed[point])-Diameter/2;
	makeOval(x,y,Diameter,Diameter);
	run("Fill", "slice");
}
// Save dots
selectImage(ID);
saveAs("tiff",pathOut+"DotsHypo");
close();


// Plot epiblast motion ////////////////////////////////////////////////////////////
// Generate blank image
newImage("HyperStack", "8-bit color-mode",ImageSize[0],ImageSize[1], 1,MaxFrame, 1);
ID=getImageID();
// Draw track after track
for (track=1;track<=NumberTracks;track++){
	count=0;
	for (i=0;i<Slicen.length;i++){
		if (Trackn[i]==track){
			// For the first tracked point, don't do anything
			count=count+1;
			if (count==1){
				X=XsEpi[i];
				Y=YsEpi[i];
			}
			if (count!=1){
				// For the 2nd to the last tracked point, draw the (n-1)-n segment
				X2=XsEpi[i];
				Y2=YsEpi[i];
				// Identify the slice and get to the slice
				slice=Slicen[i];
				// 
				for (s=slice;s<=MaxFrame;s++){
					selectImage(ID);
					setSlice(s);
					// Draw the segment
					makeLine(X,Y,X2,Y2);
					run("Draw", "slice");
				}
				// Keep the edge of the segment for drawing next segment
				X=X2;
				Y=Y2;
			}
		}
	}
}
// Save cumulated tracks
selectImage(ID);
saveAs("tiff",pathOut+"TracksEpi");
close();


// Plot epiblast motion (dots) /////////////////////////////////////////////////
// Generate blank image
newImage("HyperStack", "8-bit color-mode",ImageSize[0],ImageSize[1], 1,MaxFrame, 1);
ID=getImageID();
// Draw tracked point after tracked point
for (point=0;point<Trackn.length;point++){
	selectImage(ID);
	slice=Slicen[point];
	setSlice(slice);
	x=parseFloat(XsEpi[point])-Diameter/2;
	y=parseFloat(YsEpi[point])-Diameter/2;
	makeOval(x,y,Diameter,Diameter);
	run("Fill", "slice");
}
// Save dots
selectImage(ID);
saveAs("tiff",pathOut+"DotsEpi");
close();



// Plot differential motion over time (tracks) /////////////////////////////////////////////////
// Generate blank image
newImage("HyperStack", "8-bit color-mode",ImageSize[0],ImageSize[1], 1,MaxFrame, 1);
ID=getImageID();
// Draw track after track
for (track=1;track<=NumberTracks;track++){
	count=0;
	for (i=0;i<Slicen.length;i++){
		if (Trackn[i]==track){
			// For the first tracked point, don't do anything
			count=count+1;
			if (count==1){
				X=XsDiff[i];
				Y=YsDiff[i];
			}
			if (count!=1){
				// For the 2nd to the last tracked point, draw the (n-1)-n segment
				X2=XsDiff[i];
				Y2=YsDiff[i];
				// Identify the slice and get to the slice
				slice=Slicen[i];
				// 
				for (s=slice;s<=MaxFrame;s++){
					selectImage(ID);
					setSlice(s);
					// Draw the segment
					makeLine(X,Y,X2,Y2);
					run("Draw", "slice");
				}
				// Keep the edge of the segment for drawing next segment
				X=X2;
				Y=Y2;
			}
		}
	}
}
// Save cumulated tracks
selectImage(ID);
saveAs("tiff",pathOut+"TracksDiff");
close();


// Plot differential motion over time (dots) /////////////////////////////////////////////////
// Generate blank image
newImage("HyperStack", "8-bit color-mode",ImageSize[0],ImageSize[1], 1,MaxFrame, 1);
ID=getImageID();
// Draw tracked point after tracked point
for (point=0;point<Trackn.length;point++){
	selectImage(ID);
	slice=Slicen[point];
	setSlice(slice);
	x=parseFloat(XsDiff[point])-Diameter/2;
	y=parseFloat(YsDiff[point])-Diameter/2;
	makeOval(x,y,Diameter,Diameter);
	run("Fill", "slice");
}
// Save dots
selectImage(ID);
saveAs("tiff",pathOut+"DotsDiff");
close();
