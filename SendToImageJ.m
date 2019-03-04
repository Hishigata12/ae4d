function SendToImageJ(handles,overlay) 
    
    %Reducing the 4D image based on range
%     param = evalin('base','param');
    
    %Checks if the data set is PE or AE
    if handles.PE_4dbox.Value == 0
        if handles.hotcold.Value == 1
            if handles.bbdb.Value == 1
                h = hotcoldDB;
            else
                h = 'HotAndCold';
            end
        else
            h = 'MagnitudeColorMap';
        end
        colorMap = h;
    %    colorMap = 'MagnitudeColorMap';  
    else
        colorMap = 'grayscale';  
    end
    
    if handles.use_chop.Value == 1
        Xfilt = evalin('base','X_c');
        ax = evalin('base','ax_c');
    else
        Xfilt = evalin('base','Xfilt');
        ax = evalin('base','ax');
    end
    
    %Makes sure that the data is not complex
    Xfilt = real(Xfilt);


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

    if length(dims) == 4
        q.t = 1:dims(4);
        tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
    else
        tInd = 1;
    end
    
    %Selects the data in the specified range
    mainSet = Xfilt(xInd,yInd,zInd,tInd);

    
    %Gets the dBRng for the data
    if(isfield(handles.aeR,'String') || isempty(handles.aeR.String))  
        dBRng = [min(min(min(min(mainSet)))) max(max(max(max(mainSet))))];
    else
       dBRng = str2num(handles.aeR.String);     
    end

    
    %Checks connection of Mij
    location = checkMijConnection;
    
    sizeData = size(mainSet);
    mainSet = mainSet;
    if ndims(mainSet) == 4
    %Sends each 3D figure to imageJ
    for timePoint = 1:sizeData(4)
        %Creates the ImagePlus by changing the orientation of the data
        mainSetPartition = single(mainSet(:,:,:,timePoint));    


%         data(data<dBRng(1)) = dBRng(1);
%         data(data>dBRng(2)) = dBRng(2);

        %Creating small border so the image cannot crash even if the 3D
        %image is empty
        mainSetPartition(1,:,:) = (dBRng(2) - dBRng(1))/100 + dBRng(1);
        mainSetPartition(1,:,1) = (dBRng(2) - dBRng(1))/100 + dBRng(1);
        mainSetPartition(1,1,:) = (dBRng(2) - dBRng(1))/100 + dBRng(1);
        mainSetPartition(end,:,:) = (dBRng(2) - dBRng(1))/100 + dBRng(1);
        mainSetPartition(end,:,end) = (dBRng(2) - dBRng(1))/100 + dBRng(1);
        mainSetPartition(end,end,:) = (dBRng(2) - dBRng(1))/100 + dBRng(1);
        
        MIJ.createImage(mainSetPartition);
        size(mainSetPartition)
        
    end
    else 
        mainSetPartition = single(mainSet);
        MIJ.createImage(mainSetPartition);
    end
    
    %Writes the information that ImageJ needs to threshold and choose
    %colormap
    %Putting in fiji because it is easy to find in that directory for the
    %macro
    fileID = fopen(strcat(location,'\Fiji.app\macros\exp.txt'),'w');
    fprintf(fileID,strcat(colorMap,'\n')); %Row 0
    fprintf(fileID,strcat(num2str(dBRng(1)),'\n')); %Row 1
    fprintf(fileID,strcat(num2str(dBRng(2)),'\n')); %Row 2
    if 1==1%struc.scale.travel==1 %for PA only
        fprintf(fileID,'AE\n'); %Row 3
    else
        fprintf(fileID,'PE\n'); %Row 3
    end

    %Checks if the 4D data wants to be saved into a tif
    if handles.save_fig.Value == 1
        fprintf(fileID,'saveTrue\n'); %Row 4
    else
        fprintf(fileID,'saveFalse\n'); %Row 4
    end
    
    %Sets the name and location of the file to save
    fprintf(fileID,strcat(strjoin(strsplit(location,'\'),'\\\'), '_', '.tif','\n')); %Row 5
    
    %Simply gets the end points
    %This needs to give the correct values on the 4DCreation.ijm or 
    %4DCreationPEOverlay.ijm, needs to give the pixel width,pixel height,voxel depth
    fprintf(fileID, strcat(num2str(abs(xR(1)-xR(2))/sizeData(1)),'\n')); %Row 6
    fprintf(fileID, strcat(num2str(abs(yR(1)-yR(2))/sizeData(2)),'\n')); %Row 7
    fprintf(fileID, strcat(num2str(abs(zR(1)-zR(2))/sizeData(3)),'\n')); %Row 8
    
    %Calls the macro that turns the 3D data that was sent to ImageJ into a 4D
    %data set and puts that 4D image in the 3D viewer
    if overlay==1
        macro_path=strcat(location,'\Fiji.app\macros\4DCreationPEOverlay.ijm');
        IJObject = ij.IJ();
        IJObject.runMacroFile(macro_path);     
    else
    	macro_path=strcat(location,'\Fiji.app\macros\4DCreation.ijm');
        IJObject = ij.IJ();
        IJObject.runMacroFile(macro_path); 
    end
end