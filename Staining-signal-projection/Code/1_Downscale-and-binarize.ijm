// Downgrade and binarize image for superficial layer detection

// Parameter //////////////////////////////////////////////////////////////////
RollingBallRadius=20;
ScalingFactor=0.05;
MeanFilterSize=2;
//////////////////////////////////////////////////////////////////////////////

// CODE //////////////////////////////////////////////////////////////////////////////////////

// User select the path where stacks are stored
path = getDirectory("Select the folder containing images for topomap detection");
files = getFileList(path);

// Output folder creation
pathOut = path + "TopoMapDetection" + File.separator;
File.makeDirectory(pathOut);

// Treatment for each stack
for (i=0; i<files.length; i++) { 
		// Opening the stack
        if(endsWith(files[i],"_HOECHST.tif")) {
        	open(path+files[i]);
        	ID=getImageID();
        	selectImage(ID);
        	run("Select All");
			run("Duplicate...", "duplicate");
			IDOriginal=getImageID();
			setSlice(4);
			run("Enhance Contrast", "saturated=0.35");
			
			// Rolling ball to extract dynamic contours
			selectImage(ID);
        	run("Subtract Background...", "rolling="+RollingBallRadius+" stack");
        	
        	// Downscaling the image
        	getDimensions(width, height, channels, slices, frames);
        	newWidth=floor(width*ScalingFactor);
        	newHeight=floor(height*ScalingFactor);
        	run("Scale...", "x="+ScalingFactor+" y="+ScalingFactor+" z=1.0 width="+newWidth+" height="+newHeight+" depth="+slices+" interpolation=Bilinear average process create");
        	ID2=getImageID();
        	selectImage(ID);
        	close();
        	selectImage(ID2);
        	
        	// 2D median filter
        	run("Mean...", "radius="+MeanFilterSize+" stack");
        	run("Tile");
        	
        	// Binarize image
        	selectImage(ID2);
        	setSlice(4);
			run("Enhance Contrast", "saturated=0.35");
        	waitForUser("Apply appropriate threshold");
        	
			// Saving image
			selectImage(ID2);
			saveAs("tiff",pathOut+files[i]);
			
			// Closing image
        	selectImage(ID2);
        	close();
        	selectImage(IDOriginal);
        	close();
        	}
}