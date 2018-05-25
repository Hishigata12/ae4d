function SendToImageJ() 
    %Assigns the 4D set to the chopped set
    mainSet = evalin('base','X_c');
    
    %Checks if it is pulse echo or AE
    if 1==1 %for PA only
        dBRng    = [-12,0];
        colorMap = 'MagnitudeColorMap';  
    else
%         colorMap = 'grayscale';  
    end

    %Checks connection of Mij
    location = checkMijConnection;
    
    sizeData = size(mainSet);
    %Sends each 3D figure to imageJ
    for timePoint = 1:sizeData(4)
        %Creates the ImagePlus by changing the orientation of the data
        mainSetPartition = single(mainSet(:,:,:,timePoint));    


%         data(data<dBRng(1)) = dBRng(1);
%         data(data>dBRng(2)) = dBRng(2);

        %Creating small border so the image cannot crash even if the 3D
        %image is empty
        mainSetPartition(1,1,1) = dBRng(2);
        mainSetPartition(1,end,1) = dBRng(2);
        mainSetPartition(end,end,1) = dBRng(2);
        mainSetPartition(end,1,1) = dBRng(2);
        
        MIJ.createImage(mainSetPartition);
    end
    
    %Writes the information that ImageJ needs to threshold and choose
    %colormap
    %Putting in fiji because it is easy to find in that directory for the
    %macro
    fileID = fopen(strcat(location,'\Fiji.app\macros\exp.txt'),'w');
    fprintf(fileID,strcat(colorMap,'\n'));
    fprintf(fileID,strcat(num2str(dBRng(1)),'\n'));
    fprintf(fileID,strcat(num2str(dBRng(2)),'\n'));
    if 1==1%struc.scale.travel==1 %for PA only
        fprintf(fileID,'AE\n');
    else
        fprintf(fileID,'PE\n');
    end

    %Checks if the 4D data wants to be saved into a tif
    if 1==0
        fprintf(fileID,'saveTrue\n');
    else
        fprintf(fileID,'saveFalse\n');
    end
    
    %Sets the name and location of the file to save
    fprintf(fileID,strcat(strjoin(strsplit(location,'\'),'\\\'), '_', '.tif','\n'));
    
    %Simply gets the end points
    %This needs to give the correct values on the 4DCreation.ijm or 
    %4DCreationPEOverlay.ijm, needs to give the pixel width,pixel height,voxel depth
    fprintf(fileID, strcat(num2str(abs(-5-5)/sizeData(1)),'\n'));
    fprintf(fileID, strcat(num2str(abs(-5-5)/sizeData(2)),'\n'));
    fprintf(fileID, strcat(num2str(abs(68-80)/sizeData(3)),'\n'));
    
    %Calls the macro that turns the 3D data that was sent to ImageJ into a 4D
    %data set and puts that 4D image in the 3D viewer
    if 1==0
        macro_path=strcat(location,'\Macro\4DCreationPEOverlay.ijm');
        IJObject = ij.IJ();
        IJObject.runMacroFile(macro_path);     
    else
    	macro_path=strcat(location,'\Macro\4DCreation.ijm');
        IJObject = ij.IJ();
        IJObject.runMacroFile(macro_path); 
    end
end