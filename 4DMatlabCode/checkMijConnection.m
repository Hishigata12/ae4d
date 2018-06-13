function location = checkMijConnection()

    %Gets the location of the file by locating the location of this file
    %and back up two locations
    str = (which('checkMijConnection.m'));
    temp = strsplit(str,'\');    location = strjoin(temp(1:size(temp,2)-2),'\');
    %Checks if ImageJ is connected to MATLAB and connects if not
    if not(any(~cellfun('isempty',strfind(javaclasspath('-dynamic'),'ij.jar'))))
        fprintf('Connecting ImageJ...\n');
        %Have to download imageJ and run it to get the ij file and then put a copy
        %of that into the java folder of MATLAB
        javaaddpath(strcat(location,'\Fiji.app\ij.jar'));
        %Google mij and download mij.jar than put that into the java folder of matlab
        javaaddpath(strcat(location,'\Fiji.app\mij.jar'));
        addpath(strcat(location,'\Fiji.app\scripts'));

        %Runs MIJI in the background
        Miji(false);
        fprintf('ImageJ Connected\n');
    end
    
    %Closes all windows without saving
    while(ij.WindowManager.getCurrentImage() ~= [])
        IMG = ij.WindowManager.getCurrentImage();
        IMG.changes = false; 
        IMG.close();
    end
    MIJ.start();
end