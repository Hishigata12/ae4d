function SendToImageJ(handles) 
    
    %Reducing the 4D image based on range
    param = evalin('base','param');
    
    %Checks if the data set is PE or AE
    if handles.PE_4dbox.Value == 1
        colorMap = 'MagnitudeColorMap';  
    else
        colorMap = 'grayscale';  
    end
    %Gets the dBRng for the data
    dBRng = str2num(handles.aeR.String);
    
    if handles.use_chop.Value == 1
        Xfilt = evalin('base','X_c');
        ax = evalin('base','ax_c');
    else
        Xfilt = evalin('base','Xfilt');
        ax = evalin('base','ax');
    end


    xR = str2num(handles.xR.String);
    if length(xR) == 1
        xR = [xR xR];
    end
    yR = str2num(handles.yR.String);
    if length(yR) == 1
        yR = [yR yR];
    end
    zR = str2num(handles.zR.String);
    if length(zR) == 1
        zR = [zR zR];
    end
    tR = str2num(handles.tR.String);
    if length(tR) == 1
        tR = [tR tR];
    end
    dims = size(Xfilt);
    if length(dims) < 3
        dims(3) = 1;
    end
    %[~,ax] = make_axes(param,dims,[1 2],1);
    q.x = 1:dims(1);
    q.y = 1:dims(2);
    q.z = 1:dims(3);

    xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
    yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
    zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));

    q.t = 1:dims(4);
    tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
    
    %Selects the data in the specified range
    mainSet = Xfilt(xInd,yInd,zInd,tInd);

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
    if handles.save_fig.Value == 1
        fprintf(fileID,'saveTrue\n');
    else
        fprintf(fileID,'saveFalse\n');
    end
    
    %Sets the name and location of the file to save
    fprintf(fileID,strcat(strjoin(strsplit(location,'\'),'\\\'), '_', '.tif','\n'));
    
    %Simply gets the end points
    %This needs to give the correct values on the 4DCreation.ijm or 
    %4DCreationPEOverlay.ijm, needs to give the pixel width,pixel height,voxel depth
    fprintf(fileID, strcat(num2str(abs(xR(1)-xR(2))/sizeData(1)),'\n'));
    fprintf(fileID, strcat(num2str(abs(yR(1)-yR(2))/sizeData(2)),'\n'));
    fprintf(fileID, strcat(num2str(abs(zR(1)-zR(2))/sizeData(3)),'\n'));
    
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