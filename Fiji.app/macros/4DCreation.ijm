//Opens and splits the user created text file
filecontent=File.openAsString(getDirectory("macros")+"exp.txt");


rows=split(filecontent, "\n");
//Grabs the data from the user created text file
colorMap = rows[0];
firstBound = parseFloat(rows[1]);
secondBound = parseFloat(rows[2]);
AEorPE = rows[3];
saveTif = rows[4];
fileLocation = rows[5];
x_pixel_width = parseFloat(rows[6]);
y_pixel_width = parseFloat(rows[7]);
z_pixel_width = parseFloat(rows[8]);


//Puts the 3D images into a 4D image with the 4th Dimension being time
run("Concatenate...", "all_open title=[4D Image] open");

//Moves the window a set location
selectWindow("4D Image");
setLocation(100,100,500,500);

//Gives colormap to the 4D image
run(colorMap);

//Sets the range of the image
//run("Brightness/Contrast...");//Temp just for creator
setMinAndMax(firstBound, secondBound);


//Closes and opens 3D viewer
call("ij3d.ImageJ3DViewer.close");
run("3D Viewer");

//This will be how calibrating will be done, really useful but needs matlab to do some math
//run("Properties...", "unit=mm pixel_width=1 pixel_height=1 voxel_depth=z_pixel_width");
run("Properties...", "unit=mm pixel_width=x_pixel_width pixel_height=y_pixel_width voxel_depth=z_pixel_width");

//This is for the AE or PE part
if(AEorPE == "AE"){
    run("8-bit");
    call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false");
    call("ij3d.ImageJ3DViewer.add", "4D Image", "None", "4D Image", "0", "true", "true", "true", "2", "0");
}
else{
    run("8-bit");
    call("ij3d.ImageJ3DViewer.add", "4D Image", "None", "PE_4D_Image", "0", "true", "true", "true", "2", "0");

}

if(saveTif == "saveTrue"){
    saveAs("Tiff", fileLocation);
}

