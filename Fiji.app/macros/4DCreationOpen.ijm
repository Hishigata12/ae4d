//Opens and splits the user created text file
filecontent=File.openAsString("C:\\Users\\mynam\\Documents\\AndresMove\\Fiji.app\\macros\\exp.txt");
rows=split(filecontent, "\n");
//Grabs the data from the user created text file
colorMap = rows[0];
firstBound = parseFloat(rows[1]);
secondBound = parseFloat(rows[2]);
AEorPE = rows[3];


//run("getTransform");

//call("ij3d.ImageJ3DViewer.getUniv");
//call("ij3d.apply.close");
//print("Didn't crash new");


//This is for the AE or PE part
if(AEorPE == "AE"){
    selectWindow("4D Image");
    run(colorMap);
    run("Duplicate...", "title=[Duplicate 4D Image] duplicate");
    selectWindow("Duplicate 4D Image");
    run("8-bit");
    call("ij3d.ImageJ3DViewer.add", "Duplicate 4D Image", "None", "AE_4D_Image", "0", "true", "true", "true", "2", "0");
    close("Duplicate 4D Image");
}
else{
    call("ij3d.ImageJ3DViewer.add", "Concatenated Stacks", "None", "Concatenated_Stacks_PE", "0", "true", "true", "true", "2", "0");
    call("ij3d.ImageJ3DViewer.select", "Concatenated Stacks PE");
    call("ij3d.ImageJ3DViewer.setTransparency", "0.85");
}