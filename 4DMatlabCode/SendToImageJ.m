function SendToImageJ(struc)
    tic;
    profile clear;
    profile on;
    global struc4D;
    
    %This bit of code is here for getting the PE and AE for verasonics
    if isfield(struc,'PE') && isempty(struc.PE)~=1 && struc.p.PE
        struc.p.PE = 0;
        GettingAE4D(struc);
        global dim4PESet;
        mainSet = dim4PESet;
        struc.scale.travel = 2;
        struc.p.PE = 1;
        smoothBoolean = 1;
    else
        global dim4Set;
        smoothBoolean = 0;
        mainSet = dim4Set;
     end
    
    %Checks if it is pulse echo or AE
    if struc.scale.travel==1 %for PA only
        z_offset             = struc.scale.z_offset;
        struc.scale.dBRng    = struc.scale.dBRngCH1;
        struc4D.scale.refMAG   = struc.scale.refMAGCH1;
        colorMap = 'MagnitudeColorMap';  
    else
        z_offset             = 0;    
        struc.scale.dBRng    = struc.scale.dBRngCH2;
        struc4D.scale.refMAG   = struc.scale.refMAGCH2;
        colorMap = 'grayscale';  
    end
    
    %Doing this so 4DSet is not going to be read in everytime because it
    %is a slow process
    if(not(isfield(struc.dim4,'xROIidx')))
        %Put these here because at beginning, if started with 4D set they were undefined
        struc.dim4.xROIidx = find(struc.roi.xFile>=struc.roi.x0 & struc.roi.xFile<=struc.roi.xf);
        struc.dim4.yROIidx = find(struc.roi.yFile>=struc.roi.y0 & struc.roi.yFile<=struc.roi.yf);
        struc.dim4.zROIidx = find(struc.roi.zFile>=(struc.roi.z0-z_offset) & struc.roi.zFile<=(struc.roi.zf-z_offset));
        struc.dim4.tSlowROIidx = find(struc.roi.tSlowFile>=struc.roi.t0 & struc.roi.tSlowFile<=struc.roi.tf);
    end

    %Checks connection of Mij
    location = checkMijConnection;
    
    %Closes all windows without saving
    while(ij.WindowManager.getCurrentImage() ~= [])
        IMG = ij.WindowManager.getCurrentImage();
        IMG.changes = false; 
        IMG.close();
    end
      
    %Uses the struc when 4D data was processed to get colormap and refMAX
    %But uses current struc to get the min and max or BBRng
    if struc.scale.travel==1
        if struc4D.p.plotWFflag==1 % real

            % Check for planewaves if axial focus=0;
            if isfield(struc4D,'PE') && isempty(struc4D.PE)~=1 && struc4D.PE.PEParm.ScanParm.AxialFocus<=0   
                colorMap = 'grayscale';
                if p.BB==1
                    disp('BB Planewave reconstruction not supported (nor does t make sense) Baseband AFTER reconstrcution!');
                    return
                else
                    PWfname = [struc4D.file.saveFile '_PW_','_y' num2str(struc4D.roi.yC,2),'_S',num2str(1)];%(slice)];%Don't know what happens if I just use any slice for the range
                end
                recon            = recon_russ(PWfname);            

                dBRng = [-0.95*recon.realRef,0.95*recon.realRef];
            else
                %This one is for BaseBand and Wavefield
                colorMap = 'HotAndCold';
                bbMax = 10^(struc4D.scale.refMAG/20);
                dBRng = [-0.95*bbMax,0.95*bbMax];
            end
        elseif struc4D.p.BB==1 %This one is for BaseBand
            colorMap = 'HotAndCold2';
            dBRng = [-struc.p.dBRngCorr,struc.p.dBRngCorr];
        else
            %Check for planewaves
            dBRng = struc.scale.dBRng;
            if isfield(struc,'PE') && isempty(struc.PE)~=1 && struc.PE.PEParm.ScanParm.AxialFocus<=0 && p.BB==1
                disp('BB Planewave reconstruction not supported (nor does t make sense) Baseband AFTER reconstrcution!');
                return
            end
        end
    %Basic pulse echo dbRng
    else
        dBRng = struc.scale.dBRng;
    end

    %Needs to find the closest value to itself so if the value given is not
    %in the range of the matrix it will not return by mistake, only for the
    %checking method, if I did not need to return this and the bottom
    %checking would not be needed
    struc.roi.x0 = struc.roi.xFile(find(struc.roi.xFile>=struc.roi.x0,1,'first'));
    struc.roi.xf = struc.roi.xFile(find(struc.roi.xFile<=struc.roi.xf,1,'last'));
    struc.roi.y0 = struc.roi.yFile(find(struc.roi.yFile>=struc.roi.y0,1,'first'));
    struc.roi.yf = struc.roi.yFile(find(struc.roi.yFile<=struc.roi.yf,1,'last'));
    struc.roi.z0 = struc.roi.zFile(find(struc.roi.zFile>=struc.roi.z0,1,'first'));
    struc.roi.zf = struc.roi.zFile(find(struc.roi.zFile<=struc.roi.zf,1,'last'));
    struc.roi.t0 = struc.roi.tSlowFile(find(struc.roi.tSlowFile>=struc.roi.t0,1,'first'));
    struc.roi.tf = struc.roi.tSlowFile(find(struc.roi.tSlowFile<=struc.roi.tf,1,'last'));
    
    
    %Gets the new range needed
    %If the user gives an out of bound range the system quits and tells the
    %user the part that is out of bounds
    sizeData = size(mainSet);
    xROIidx = struc.roi.xFile(struc.dim4.xROIidx);
    if(xROIidx(1)>struc.roi.x0 || xROIidx(end)<struc.roi.xf)
        fprintf('ROI is out of range for the previous 4D set render on X-axis\n');
        fprintf('4D Range is from %f to %f\n',xROIidx(1),xROIidx(end));
        return;
    end
    interpIndex = linspace(xROIidx(1),xROIidx(end),sizeData(2));
    xROIidx = find(interpIndex>=struc.roi.x0 & interpIndex<=struc.roi.xf);
    yROIidx = struc.roi.yFile(struc.dim4.yROIidx);
    if(yROIidx(1)>struc.roi.y0 || yROIidx(end)<struc.roi.yf)
        fprintf('ROI is out of range for the previous 4D set render on Y-axis\n');
        fprintf('4D Range is from %f to %f\n',yROIidx(1),yROIidx(end));
        return;
    end
    interpIndex = linspace(yROIidx(1),yROIidx(end),sizeData(1));
    yROIidx = find(interpIndex>=struc.roi.y0 & interpIndex<=struc.roi.yf);
    zROIidx = struc.roi.zFile(struc.dim4.zROIidx);
    if(zROIidx(1)>struc.roi.z0 || zROIidx(end)<struc.roi.zf)
        fprintf('ROI is out of range for the previous 4D set render on Z-axis\n');
        fprintf('4D Range is from %f to %f\n',zROIidx(1),zROIidx(end));
        return;
    end
    interpIndex = linspace(zROIidx(1),zROIidx(end),sizeData(3));
    zROIidx     = find(interpIndex>=(struc.roi.z0-struc.scale.z_offset) & interpIndex<=(struc.roi.zf-struc.scale.z_offset));
    tSlowROIidx = struc.roi.tSlowFile(struc.dim4.tSlowROIidx);
    if(tSlowROIidx(1)>struc.roi.t0 || tSlowROIidx(end)<struc.roi.tf)
        fprintf('ROI is out of range for the previous 4D set render on T-axis\n');
        fprintf('4D Range is from %f to %f\n',tSlowROIidx(1),tSlowROIidx(end));
        return;
    end
    interpIndex = linspace(tSlowROIidx(1),tSlowROIidx(end),sizeData(4));
    tSlowROIidx = find(interpIndex>=struc.roi.t0 & interpIndex<=struc.roi.tf);
    
    %Width is y direction, height is x direction
    %Sends each 3D figure to imageJ
    %minVal = min(mainSet(:));
    %maxVal = max(mainSet(:));
    for timePoint = tSlowROIidx
        %Creates the ImagePlus by changing the orientation of the data
        data = single(mainSet(yROIidx, xROIidx,zROIidx,timePoint));    


        data(data<dBRng(1)) = dBRng(1);
        data(data>dBRng(2)) = dBRng(2);
        
        %Smoothing data
        w=7;
        smoothZData = data;
        if smoothBoolean==1
            sum = zeros(size(data,1),size(data,2));
            for k = (w+1):size(data,3)-w
                for j = (-w+k):k+w
                    sum = sum + squeeze(data(:,:,j));
                end
                smoothZData(:,:,k) = sum/(2*w+1);
                sum(:,:) = 0;
            end
        end
        %Creating small border so the image cannot crash even if the 3D
        %image is empty
        smoothZData(1,1,1) = dBRng(2);
        smoothZData(1,end,1) = dBRng(2);
        smoothZData(end,end,1) = dBRng(2);
        smoothZData(end,1,1) = dBRng(2);
        
        MIJ.createImage(smoothZData);
    end
    
    %Writes the information that ImageJ needs to threshold and choose
    %colormap
    %Putting in fiji because it is easy to find in that directory for the
    %macro
    fileID = fopen(strcat(location,'\Fiji.app\macros\exp.txt'),'w');
    fprintf(fileID,strcat(colorMap,'\n'));
    fprintf(fileID,strcat(num2str(dBRng(1)),'\n'));
    fprintf(fileID,strcat(num2str(dBRng(2)),'\n'));
    if struc.scale.travel==1 %for PA only
        fprintf(fileID,'AE\n');
    else
        fprintf(fileID,'PE\n');
    end

    %Checks if the 4D data wants to be saved into a tif
    if struc.file.saveData3D
        fprintf(fileID,'saveTrue\n');
    else
        fprintf(fileID,'saveFalse\n');
    end
    
    %Sets the name and location of the file to save
    fprintf(fileID,strcat(strjoin(strsplit(struc.file.saveFile,'\'),'\\\'), '_', '.tif','\n'));
    
    %Simply gets the end points
    %This needs to give the correct values on the 4DCreation.ijm or 
    %4DCreationPEOverlay.ijm, needs to give the pixel width,pixel height,voxel depth
    fprintf(fileID, strcat(num2str(abs((struc.roi.x0-struc.roi.xf))/sizeData(2)),'\n'));
    fprintf(fileID, strcat(num2str(abs((struc.roi.y0-struc.roi.yf))/sizeData(1)),'\n'));
    fprintf(fileID, strcat(num2str(abs((struc.roi.z0-struc.roi.zf))/sizeData(3)),'\n'));
    
    %Calls the macro that turns the 3D data that was sent to ImageJ into a 4D
    %data set and puts that 4D image in the 3D viewer
    if isfield(struc,'PE') && isempty(struc.PE)~=1 && struc.p.PE
        macro_path=strcat(location,'\macros\4DCreationPEOverlay.ijm');
        IJObject = ij.IJ();
        IJObject.runMacroFile(macro_path);     
    else
    	macro_path=strcat(location,'\macros\4DCreation.ijm');
        IJObject = ij.IJ();
        IJObject.runMacroFile(macro_path); 
    end
    
% macro_path=strcat(location,'\macros\4DCreation.ijm');
% IJObject = ij.IJ();
% IJObject.runMacroFile(macro_path); 
%         
profile off;
profile viewer;
toc;
end