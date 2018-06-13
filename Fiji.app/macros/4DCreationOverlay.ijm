//Puts the 3D images into a 4D image with the 4th Dimension being time
run("Concatenate...", "all_open title=[4D Image] open");

//Opens and splits the user created text file
filecontent=File.openAsString("C:\\Users\\mynam\\Documents\\AndresMove\\Fiji.app\\macros\\exp.txt");
rows=split(filecontent, "\n");
//Grabs the data from the user created text file
colorMap = rows[0];
firstBound = parseFloat(rows[1]);
secondBound = parseFloat(rows[2]);
AEorPE = rows[3];

//Moves the window a set location
selectWindow("4D Image");
setLocation(500,500,500,500);

//Gives colormap to the 4D image
run(colorMap);

//Sets the range of the image
//run("Brightness/Contrast...");//Temp just for creator
setMinAndMax(firstBound, secondBound);