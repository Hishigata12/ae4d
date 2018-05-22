//Opens and splits the user created text file
//close();
filecontent=File.openAsString(getDirectory("macros")+"exp.txt");
rows=split(filecontent, "\n");
//Grabs the data from the user created text file
AETifLocation = rows[0];
PETifLocation = rows[1];

//Open and closes the image3D window
call("ij3d.ImageJ3DViewer.close");
run("3D Viewer");

//Opens the AE Tif image and gets dimensions
open(AETifLocation);
splitFile =  split(AETifLocation,"\\");
call("ij3d.ImageJ3DViewer.add", splitFile[splitFile.length-1], "None", "AE_4D_Image", "0", "true", "true", "true", "2", "0");
getDimensions(wAE, hAE, channelsAE, slicesAE, framesAE);

//Opens the PE Tif image, gets dimensions and puts a transparency
open(PETifLocation);
splitFile =  split(PETifLocation,"\\");
getDimensions(wPE, hPE, channelsPE, slicesPE, framesPE);
run("Size...", "width=wAE height=hAE depth=slicesAE time=framesAE average interpolation=Bilinear");
run("8-bit");
call("ij3d.ImageJ3DViewer.add", splitFile[splitFile.length-1], "None", "PE_4D_Image", "0", "true", "true", "true", "2", "0");
call("ij3d.ImageJ3DViewer.select", "PE_4D_Image");
call("ij3d.ImageJ3DViewer.setTransparency", "0.1");

//close();
