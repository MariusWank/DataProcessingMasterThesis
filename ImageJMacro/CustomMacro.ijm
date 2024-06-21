// ImageJ Macro Script to analyze fractal dimensionality of multiple images and average the results

// Folder containing images
inputDir = getDirectory("Choose Directory with Images");

// List of image files in the directory
fileList = getFileList(inputDir);

// Initialize variables to store total fractal dimension and count
totalFractalDimension = 0;
imageCount = 0;
runCount = 0;

// Loop through all run folders within the parameter set folder
for (runIndex = 0; runIndex < fileList.length; runIndex++) {
    // Get the current run folder
    runFolder = inputDir + fileList[runIndex];

    // Skip the folder named "measurements"
    if (fileList[runIndex] == "measurements/") {
        continue;
    }
    
    // List of PNG files within the run folder
    pngFiles = getFileList(runFolder);

    // Loop through all PNG files within the run folder
    for (pngIndex = 0; pngIndex < pngFiles.length; pngIndex++) {

        // Check if the PNG file starts with "000600"
        if (pngFiles[pngIndex].startsWith("000600")) {
            // Open the PNG image
            open(runFolder + "/" + pngFiles[pngIndex]);

            // Ensure the image is binary (convert if necessary)
            run("8-bit");
            setAutoThreshold("Default");
            run("Convert to Mask");

            // Calculate fractal dimension using box-counting method
            run("Fractal Box Count...", "box=1,2,3,4,6,8,12,16,32,64,128,256,512,1024 black");

            // Get the fractal dimension from the Results table
            fractalDimension = getResult("D", pngIndex-1);
        
            // Add the fractal dimension to the total
            totalFractalDimension += fractalDimension;
            imageCount++;

            // Close the current image
            close();

            selectImage("000600k000000.png");
            close();
        }
    }
}

// Calculate the average fractal dimension
if (imageCount > 0) {
    averageFractalDimension = totalFractalDimension / imageCount;
    print("Average Fractal Dimension: " + averageFractalDimension);
    saveAs("Results", inputDir + "/fractal_dimension.txt");
} else {
    print("No images processed.");
}
// Close all open images
while (nImages > 0) {
    selectImage(nImages);
    close();
}
