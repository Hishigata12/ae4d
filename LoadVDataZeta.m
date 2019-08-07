function [HF, LF, PE, Parms, LF1, Location] = LoadVDataZeta(Parms,ScanPt,XPt,YPt,chIdx,DataToLoad,Average,LoadPE,LF1,Location,AvgExist,UnAvgExist,UnShiftExist)
%% LoadVDataBeta.m
%  Loads high and low frequency voltage data from a scan taken using BScan.m
%  Written by Alex Alvarez 12/08/18.
%  Modified latest on 04/19/19
%  Arguments
%   Parms (Struc) - Structure of Parameters from Scan (usually output in BScan.m and then used for AnalyzeAEI.m
%   DataToLoad (String) - Will signal to load 'UnShift,' 'UnAvg,' 'Avg,' or 'Filt' data
%   Average (Boolean) - Flag to determine if the user wants to Average (or ReAverage) any Post Averaged Data
%   LoadPE (Boolean) - Flag to determine if the user wants to load PE Data as well
%   LF1 (Vector of Doubles or []) - Template Waveform for Shifting
%   Location (double or []) - Quarter Point of Location
%   AvgExist (Boolean) - Flag for if '_Raw4dAveraged.h5' file exists
%   UnAvgExist (Boolean) - Flag for if '_Raw4dUnAveraged.h5' file exists
%   UnShiftExist (Boolean) - Flag for if '_Raw4dUnShifted.h5' file exists
%  Outputs
%   HF (2-3D array of doubles) - HF data(Fast Time x Slow Time x Average)
%   LF (1-2D array of doubles) - LF data (Slow Time x Average)
%   PE (3-4D array of doubles) - PE_RF data (Fast Time x Slow Time x Element x Average)
%   Parms (struc) - Array of all parameters used in Scanning
%   LF1 (Vector of Doubles or []) - Template Waveform for Shifting
%   Location (double or []) - Quarter Point of Location
%% TOC
%   @00 High Level Flags
%   @01 Initial Definitions
%   @02 Read HF/LF Data and Average
%% @00 High Level Flags
% High level flags to control functionality of script for easier user use
LoadTimePerPt = tic;
FilterCuts = [10 30 990 1000]; % 4x1 that defines the Hanning filter cutoffs for filtering before cross-correlation
Manually = 0; % Flag to control if the user wants to manually select averages to throw away or split
ShowPlots = 0; % Flag to control if PostAverage plots (pre-shifted, correlations, post-shifted, and averaged) are shown
%% @01 Initital Definitions
% Initial definitions of quantities to be used in this script and saved to
% the Parameters files
Parms.Chan.HFChannelNumbers = str2num(strrep(Parms.Daq.HF.Channels,':',' ')); % Determine the HF Channel Numbers
% If the total number of channels is greater than 1, you will read in the
% array and convert it into a numerical array
if numel(Parms.Chan.HFChannelNumbers) > 1
    Parms.Chan.HFChanns = Parms.Chan.HFChannelNumbers(1):1:Parms.Chan.HFChannelNumbers(2); Parms.Chan.NumHFChanns = numel(Parms.Chan.HFChanns); % Set the number and names of the HF Channels
else
    Parms.Chan.HFChanns = Parms.Chan.HFChannelNumbers; Parms.Chan.NumHFChanns = 1;
end
% Repeat for the LF Channels
Parms.Chan.LFChannelNumbers = str2num(strrep(Parms.Daq.LF.Channels,':',' ')); % Determine the LF channel numbers used
if numel(Parms.Chan.LFChannelNumbers) > 1
    Parms.Chan.LFChanns = Parms.Chan.LFChannelNumbers(1):1:Parms.Chan.LFChannelNumbers(2); Parms.Chan.NumLFChanns = numel(Parms.Chan.LFChanns); % Set the number and names of the LF Channels
else
    Parms.Chan.LFChanns = Parms.Chan.LFChannelNumbers; Parms.Chan.NumLFChanns = 1;
end
NumHFChanns = Parms.Chan.NumHFChanns; % Set the number and names of the HF Channels
NumLFChanns = Parms.Chan.NumLFChanns; % Set the number and names of the LF Channels

%%
chanNum = 1:NumLFChanns; % Define an array of channels from 1 to the number of LF channels to iterate over 
Xsteps = Parms.Scan.Xpt;
Ysteps = Parms.Scan.Ypt;
if LoadPE
    cd ..\PEData % change to the appropriate path location for PE data
    DisStruc = load([Parms.Format.filename2 '_PEParm.mat']); % Load the path
    PEPRF = DisStruc.bScanParm.MaxPRF;
    Receive = DisStruc.Rcv;
    Trans = DisStruc.Trans;
    AcqsPerFrame = DisStruc.bScanParm.NumyFrames; % [double] Number of acquisition pulses per frame
    if DisStruc.bScanParm.twoAcq
        AcqsPerFrame = floor(AcqsPerFrame/2);
    end
    Diffie = Receive(2).Apod+Receive(1).Apod;
    if max(Diffie(:)) > 1
        NumElz = sum(Receive(1).Apod);
    else
        NumElz = sum(Diffie);
    end
    AcqSize = Receive(1).endSample - Receive(1).startSample + 1; % [double] Size of each pulse acquisition
% %     DataSize = size(RF_Data); % [doubles] Vector of sizes for PEMMode     
    cd ..\ExpData
    clear DisStruc;
else
    PEPRF = [];
    NumElz = 1;
end
idx = ScanPt; % (YPt-1)*Xsteps+XPt; % Define idx as you did with the scan point above
tmp1 = [pwd, '\', Parms.Format.filename2, '_Raw4dUnShifted.h5']; % Create a filename for saving (and a directory)
tmp2 = '/LF'; % Create a temporary variable for an h5 saving variable
tmp3 = '/HF'; % Repeat with HF
tmp4 = [pwd, '\', Parms.Format.filename2, '_Raw4dUnAveraged.h5'];
tmp5 = [pwd, '\', Parms.Format.filename2, '_Raw4dAveraged.h5'];
tmp6 = '/PE'; % Repeat with PE
tmp7 = [pwd, '\', Parms.Format.filename2, '_Raw4d.h5'];

sLFUnAvg =[Parms.Daq.LF.Samples,Parms.Scan.Avg,NumLFChanns,Inf,Inf]; % Set the size of the LF array that you will save for the UnAveraged Data
csLFUnAvg = [ceil(sLFUnAvg(1)/10),ceil(sLFUnAvg(2)/10),ceil(sLFUnAvg(3)/10),ceil(Xsteps/10),ceil(Ysteps/10)]; % Set the chunk size (for memory use)
sHFUnAvg =[round(Parms.daq.HFdaq.pts),Parms.daq.HFdaq.NoBurstTriggers,Parms.Scan.Avg,NumHFChanns,Inf,Inf];  % Set the size of the LF array that you will save for the UnAveraged Data
csHFUnAvg = [ceil(sHFUnAvg(1)/10),ceil(sHFUnAvg(2)/10),ceil(sHFUnAvg(3)/10),ceil(sHFUnAvg(4)/10),ceil(Xsteps/10),ceil(Ysteps/10)]; % Set the chunk size (for memory use)
sPEUnAvg =[AcqSize, AcqsPerFrame,NumElz,Parms.Scan.Avg, Inf, Inf];  % Set the size of the LF array that you will save for the UnAveraged Data
csPEUnAvg = [ceil(sPEUnAvg(1)/10),ceil(sPEUnAvg(2)/10),ceil(sPEUnAvg(3)/10),ceil(sPEUnAvg(4)/10),ceil(Xsteps/10),ceil(Ysteps/10)]; % Set the chunk size (for memory use)
sLFAvg =[Parms.Daq.LF.Samples,NumLFChanns,Inf,Inf]; % Set the size of the LF array that you will save for the UnAveraged Data
csLFAvg = [ceil(sLFAvg(1)/10), ceil(sLFAvg(2)/10),ceil(Xsteps/10),ceil(Ysteps/10)]; % Set the chunk size (for memory use)
sHFAvg =[round(Parms.daq.HFdaq.pts),Parms.daq.HFdaq.NoBurstTriggers,NumHFChanns,Inf,Inf];  % Set the size of the LF array that you will save for the UnAveraged Data
csHFAvg = [ceil(sHFAvg(1)/10),ceil(sHFAvg(2)/10),ceil(sHFAvg(3)/10),ceil(Xsteps/10),ceil(Ysteps/10)]; % Set the chunk size (for memory use)
sPEAvg =[AcqSize,AcqsPerFrame,NumElz,Inf,Inf];  % Set the size of the LF array that you will save for the UnAveraged Data
csPEAvg = [ceil(sPEAvg(1)/10),ceil(sPEAvg(2)/10),ceil(sPEAvg(3)/10),ceil(Xsteps/10),ceil(Ysteps/10)]; % Set the chunk size (for memory use)
startLFUnShift = [1 1 1 XPt YPt];
startHFUnShift = [1 1 1 1 XPt YPt];
startLFUnAvg = [1 1 chIdx XPt YPt];
startHFUnAvg = [1 1 1 chIdx XPt YPt];
startPEUnAvg  = [1 1 1 1 XPt YPt];
startLFAvg = [1 chIdx XPt YPt];
startHFAvg = [1 1 chIdx XPt YPt];
startPEAvg  = [1 1 1 XPt YPt];
countLFUnShift = [sLFUnAvg(1),sLFUnAvg(2),sLFUnAvg(3),1,1];
countHFUnShift = [sHFUnAvg(1),sHFUnAvg(2),sHFUnAvg(3),sHFUnAvg(4),1,1];
countLFUnAvg = [sLFUnAvg(1),sLFUnAvg(2),1,1,1];
countHFUnAvg = [sHFUnAvg(1),sHFUnAvg(2),sHFUnAvg(3),1,1,1];
countPEUnAvg = [sPEUnAvg(1),sPEUnAvg(2),sPEUnAvg(3),sPEUnAvg(4),1,1];
countLFAvg = [sLFAvg(1),1,1,1];
countHFAvg = [sHFAvg(1),sHFAvg(2),1,1,1];
countPEAvg = [sPEAvg(1),sPEAvg(2),sPEAvg(3),1,1];      

if strcmpi(DataToLoad,'Avg')        
    if ~isfield(Parms.Format,'Save_Avg') || ~Parms.Format.Save_Avg % Check if this is a scan where save average was checked
        if AvgExist
            if ~Average
                LF = h5read(tmp5,tmp2,startLFAvg,countLFAvg);
                HF = h5read(tmp5,tmp3,startHFAvg,countHFAvg);
                if LoadPE
                    PE = h5read(tmp5,tmp6,startPEAvg,countPEAvg);
                end
                LF1 = []; Location = [];
                return
            else
                if chIdx == 1 && ScanPt == 1
                    delete(tmp5); delete(tmp4);
                    h5create(tmp4,tmp2,sLFUnAvg,'Datatype','double','ChunkSize', csLFUnAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                    h5create(tmp4,tmp3,sHFUnAvg,'Datatype','double','ChunkSize',csHFUnAvg,'Deflate',5);
                    h5create(tmp5,tmp2,sLFAvg,'Datatype','double','ChunkSize',csLFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                    h5create(tmp5,tmp3,sHFAvg,'Datatype','double','ChunkSize',csHFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)                   
                    if LoadPE
                        h5create(tmp4,tmp6,sPEUnAvg,'Datatype','double','ChunkSize',csPEUnAvg,'Deflate',5);
                        h5create(tmp5,tmp6,sPEAvg,'Datatype','double','ChunkSize',csPEAvg,'Deflate',5); 
                    end  
                end
                LFUnShift = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
                HFUnShift = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
                if LoadPE
                    PEUnShift = h5read(tmp4,tmp6,startPEUnAvg,countPEUnAvg);
                end
            end
        elseif UnAvgExist
            if ~Average
                if chIdx == 1 && ScanPt == 1
                    h5create(tmp5,tmp2,sLFAvg,'Datatype','double','ChunkSize',csLFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                    h5create(tmp5,tmp3,sHFAvg,'Datatype','double','ChunkSize',csHFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)                   
                    if LoadPE
                        h5create(tmp5,tmp6,sPEAvg,'Datatype','double','ChunkSize',csPEAvg,'Deflate',5); 
                    end
                end
                LFUnAvg = h5read(tmp4,tmp2,startLFUnAvg,countLFUnAvg);
                HFUnAvg = h5read(tmp4,tmp3,startHFUnAvg,countHFUnAvg);
                if LoadPE
                    PEUnAvg = h5read(tmp4,tmp6,startPEUnAvg,countPEUnAvg);
                end
                LF = mean(LFUnAvg,2);
                HF = mean(HFUnAvg,3);
                PE = mean(PEUnAvg,4);
                LF1 = []; Location = [];
                h5write(tmp5,tmp2,LF,startLFAvg,countLFAvg)
                h5write(tmp5,tmp3,HF,startHFAvg,countHFAvg)
                if chIdx == 1
                    h5write(tmp5,tmp6,PE,startPEAvg,countPEAvg)
                end
                return
            else
                if chIdx == 1 && ScanPt == 1
                    delete(tmp4);
                    h5create(tmp4,tmp2,sLFUnAvg,'Datatype','double','ChunkSize', csLFUnAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                    h5create(tmp4,tmp3,sHFUnAvg,'Datatype','double','ChunkSize',csHFUnAvg,'Deflate',5);
                    h5create(tmp5,tmp2,sLFAvg,'Datatype','double','ChunkSize',csLFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                    h5create(tmp5,tmp3,sHFAvg,'Datatype','double','ChunkSize',csHFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)                
                    if LoadPE
                        h5create(tmp4,tmp6,sPEUnAvg,'Datatype','double','ChunkSize',csPEUnAvg,'Deflate',5);
                        h5create(tmp5,tmp6,sPEAvg,'Datatype','double','ChunkSize',csPEAvg,'Deflate',5); 
                    end
                end
                LFUnShift = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
                HFUnShift = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
                if LoadPE
                    PEUnShift = h5read(tmp1,tmp6,startPEUnAvg,countPEUnAvg);
                end
            end
        elseif UnShiftExist
            if chIdx == 1 && ScanPt == 1    
                h5create(tmp4,tmp2,sLFUnAvg,'Datatype','double','ChunkSize', csLFUnAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                h5create(tmp4,tmp3,sHFUnAvg,'Datatype','double','ChunkSize',csHFUnAvg,'Deflate',5);
                h5create(tmp5,tmp2,sLFAvg,'Datatype','double','ChunkSize',csLFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                h5create(tmp5,tmp3,sHFAvg,'Datatype','double','ChunkSize',csHFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                if LoadPE
                    h5create(tmp4,tmp6,sPEUnAvg,'Datatype','double','ChunkSize',csPEUnAvg,'Deflate',5);
                    h5create(tmp5,tmp6,[sPEUnAvg(1),sPEUnAvg(2),sPEUnAvg(3),sPEUnAvg(5),sPEUnAvg(6)],'Datatype','double','ChunkSize',[csPEUnAvg(1),csPEUnAvg(2),csPEUnAvg(3),csPEUnAvg(5),csPEUnAvg(6)],'Deflate',5); 
                end
            end
            LFUnShift = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
            HFUnShift = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
            if LoadPE
                PEUnShift = h5read(tmp1,tmp6,startPEUnAvg,countPEUnAvg);
            end                
        else
            if chIdx == 1
                if ScanPt == 1
                    h5create(tmp1,tmp2,sLFUnAvg,'Datatype','double','ChunkSize', csLFUnAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                    h5create(tmp1,tmp3,sHFUnAvg,'Datatype','double','ChunkSize',csHFUnAvg,'Deflate',5);
                    h5create(tmp4,tmp2,sLFUnAvg,'Datatype','double','ChunkSize', csLFUnAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                    h5create(tmp4,tmp3,sHFUnAvg,'Datatype','double','ChunkSize',csHFUnAvg,'Deflate',5);
                    h5create(tmp5,tmp2,sLFAvg,'Datatype','double','ChunkSize',csLFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                    h5create(tmp5,tmp3,sHFAvg,'Datatype','double','ChunkSize',csHFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                    % [Fast Time,SlowTime (i.e., Number of Pulses), US Array Elements, Averages, XPt, YPt]
                    if LoadPE
                        h5create(tmp1,tmp6,sPEUnAvg,'Datatype','double','ChunkSize',csPEUnAvg,'Deflate',5);
                        h5create(tmp4,tmp6,sPEUnAvg,'Datatype','double','ChunkSize',csPEUnAvg,'Deflate',5);
                        h5create(tmp5,tmp6,sPEAvg,'Datatype','double','ChunkSize',csPEAvg,'Deflate',5);
                    end    
                end
                [HFUnShift,LFUnShift,PEUnShift] = UnAverageReadZeta(Parms,idx,XPt,YPt,LoadPE); % Use the UnAverage reader if the scans have not been averaged           
                SaveUnAvgTime = tic;
                h5write(tmp1,tmp2,LFUnShift,startLFUnShift,countLFUnShift); % Write the LF variable to the h5 file
                h5write(tmp1,tmp3,HFUnShift,startHFUnShift,countHFUnShift);
                if LoadPE
                    h5write(tmp1,tmp6,PEUnShift,startPEUnAvg,countPEUnAvg); 
                end
                disp([num2str(toc(SaveUnAvgTime)),' Seconds to Save UnShifted Scan @ X=' num2str(XPt) ' Y=' num2str(YPt)])            
                LFUnShift = squeeze(LFUnShift(:,:,chIdx));
                HFUnShift = squeeze(HFUnShift(:,:,:,chIdx));
            else
                LFUnShift = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
                HFUnShift = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
                if LoadPE
                    PEUnShift = h5read(tmp1,tmp6,startPEUnAvg,countPEUnAvg);
                end                           
            end
        end
        if max(HFUnShift(:)) ~= 0 && max(LFUnShift(:)) ~=0
            % Use PostAvg to cross-correlate, shift, and average
            % both HF and LF data sets
            [HF,LF,PE,LF1,Location,NumThrown,NumSplit] = PostAverageZeta_v2(Parms,LF1,Location,LFUnShift,HFUnShift,FilterCuts,idx,XPt,YPt,chanNum(chIdx),chIdx,Manually,ShowPlots,DataToLoad,LoadPE,PEPRF,PEUnShift);                       
            % Clear the HF and LF variables from above
        else
            HF = zeros(size(HFUnShift,1),size(HFUnShift,2));
            LF = zeros(size(LFUnShift,1),1);
            if LoadPE && ~isempty(PEUnShift) && chIdx == 1
                PE = zeros(size(PEUnShift,1),size(PEUnShift,2),size(PEUnShift,3));
            end
            NumThrown = 0; NumSplit = 0;
        end
        % Determine the number of averages thrown and split in
        % each scan
        AvgNumThrown(idx) = NumThrown;
        AvgNumSplit(idx) = NumSplit;            
        SaveAvgTime = tic;
        h5write(tmp5,tmp2,LF,startLFAvg,countLFAvg); % Write the LF variable to the h5 file
        h5write(tmp5,tmp3,HF,startHFAvg,countHFAvg);
        disp([num2str(toc(SaveAvgTime)),' Seconds to Save Shifted and Averaged Data for chIdx=' num2str(chIdx) ' @ X=' num2str(XPt) ' Y=' num2str(YPt)])                    
        if LoadPE  && chIdx == 1 
            h5write(tmp5,tmp6,PE,startPEAvg,countPEAvg);
        end
    else
        if chIdx == 1 && ScanPt == 1
            h5create(tmp7,tmp2,sLFAvg,'Datatype','double','ChunkSize',csLFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
            h5create(tmp7,tmp3,sHFAvg,'Datatype','double','ChunkSize',csHFAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
            % [Fast Time,SlowTime (i.e., Number of Pulses), US Array Elements, Averages, XPt, YPt]
            if LoadPE
               h5create(tmp7,tmp6,sPEAvg,'Datatype','double','ChunkSize',csPEAvg,'Deflate',5);
            end    
        end
        [HF,LF,PE] = AverageReadZeta(Parms,idx,XPt,YPt,LoadPE);
        LF = squeeze(LF(:,chIdx));
        HF = squeeze(HF(:,:,chIdx));
        LF1 = []; Location = [];
        SaveUnShiftTime = tic;
        h5write(tmp7,tmp2,LF,startLFAvg,countLFAvg); % Write the LF variable to the h5 file
        h5write(tmp7,tmp3,HF,startHFAvg,countHFAvg);
        if LoadPE && chIdx == 1
            h5write(tmp7,tmp6,PE,startPEAvg,countPEAvg); 
        end
        disp([num2str(toc(SaveUnShiftTime)),' Seconds to Save UnShifted Data for chIdx=' num2str(chIdx) ' @ X=' num2str(XPt) ' Y=' num2str(YPt)])            
    end
elseif strcmpi(DataToLoad,'UnAvg')
    if AvgExist
        delete(tmp5)
        LF = h5read(tmp4,tmp2,startLFUnAvg,countLFUnAvg);
        HF = h5read(tmp4,tmp3,startHFUnAvg,countHFUnAvg);
        if LoadPE
            PE = h5read(tmp4,tmp6,startPEUnAvg,countPEUnAvg);
        end
        LF1 = []; Location = [];
        return        
    elseif UnAvgExist
        LF = h5read(tmp4,tmp2,startLFUnAvg,countLFUnAvg);
        HF = h5read(tmp4,tmp3,startHFUnAvg,countHFUnAvg);
        if LoadPE
            PE = h5read(tmp4,tmp6,startPEUnAvg,countPEUnAvg);
        end
        LF1 = []; Location = [];
        return
    elseif UnShiftExist
        if chIdx == 1 && ScanPt == 1    
            h5create(tmp4,tmp2,sLFUnAvg,'Datatype','double','ChunkSize', csLFUnAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
            h5create(tmp4,tmp3,sHFUnAvg,'Datatype','double','ChunkSize',csHFUnAvg,'Deflate',5);
            if LoadPE
                h5create(tmp4,tmp6,sPEUnAvg,'Datatype','double','ChunkSize',csPEUnAvg,'Deflate',5);
            end
        end
        LFUnShift = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
        HFUnShift = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
        if LoadPE
            PEUnShift = h5read(tmp1,tmp6,startPEUnAvg,countPEUnAvg);
        end                
    else
        if chIdx == 1
            if ScanPt == 1
                h5create(tmp1,tmp2,sLFUnAvg,'Datatype','double','ChunkSize', csLFUnAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                h5create(tmp1,tmp3,sHFUnAvg,'Datatype','double','ChunkSize',csHFUnAvg,'Deflate',5);
                h5create(tmp4,tmp2,sLFUnAvg,'Datatype','double','ChunkSize', csLFUnAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                h5create(tmp4,tmp3,sHFUnAvg,'Datatype','double','ChunkSize',csHFUnAvg,'Deflate',5);
                % [Fast Time,SlowTime (i.e., Number of Pulses), US Array Elements, Averages, XPt, YPt]
                if LoadPE
                    h5create(tmp1,tmp6,sPEUnAvg,'Datatype','double','ChunkSize',csPEUnAvg,'Deflate',5);
                    h5create(tmp4,tmp6,sPEUnAvg,'Datatype','double','ChunkSize',csPEUnAvg,'Deflate',5);
                end    
            end
            [HFUnShift,LFUnShift,PEUnShift] = UnAverageReadZeta(Parms,idx,XPt,YPt,LoadPE); % Use the UnAverage reader if the scans have not been averaged           
            tic
            h5write(tmp1,tmp2,LFUnShift,startLFUnShift,countLFUnShift); % Write the LF variable to the h5 file
            h5write(tmp1,tmp3,HFUnShift,startHFUnShift,countHFUnShift);
            if LoadPE
                h5write(tmp1,tmp6,PEUnShift,startPEUnAvg,countPEUnAvg); 
            end
            disp([num2str(toc),' Seconds to Save UnShifted Scan @ X=' num2str(XPt) ' Y=' num2str(YPt)])            
            LFUnShift = squeeze(LFUnShift(:,:,chIdx));
            HFUnShift = squeeze(HFUnShift(:,:,:,chIdx));
        else
            LFUnShift = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
            HFUnShift = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
            if LoadPE
                PEUnShift = h5read(tmp1,tmp6,startPEUnAvg,countPEUnAvg);
            end                           
        end
    end
    if max(HFUnShift(:)) ~= 0 && max(LFUnShift(:)) ~=0
        % Use PostAvg to cross-correlate, shift, and average
        % both HF and LF data sets
        [HF,LF,PE,LF1,Location,NumThrown,NumSplit] = PostAverageZeta_v2(Parms,LF1,Location,LFUnShift,HFUnShift,FilterCuts,idx,XPt,YPt,chanNum(chIdx),chIdx,Manually,ShowPlots,DataToLoad,LoadPE,PEPRF,PEUnShift);                       
        % Clear the HF and LF variables from above
    else
        HF = zeros(size(HFUnShift,1),size(HFUnShift,2));
        LF = zeros(size(LFUnShift,1),1);
        if LoadPE && ~isempty(PEUnShift) && chIdx == 1
            PE = zeros(size(PEUnShift,1),size(PEUnShift,2),size(PEUnShift,3));
        end
        NumThrown = 0; NumSplit = 0;
    end
    % Determine the number of averages thrown and split in
    % each scan
    AvgNumThrown(idx) = NumThrown;
    AvgNumSplit(idx) = NumSplit;            
    tic
    h5write(tmp4,tmp2,LF,startLFUnAvg,countLFUnAvg); % Write the LF variable to the h5 file
    h5write(tmp4,tmp3,HF,startHFUnAvg,countHFUnAvg);
    disp([num2str(toc),' Seconds to Save Shifted and Averaged Scan for Chan=' num2str(chanNum(chIdx)) ' @ X=' num2str(XPt) ' Y=' num2str(YPt)])                    
    if LoadPE  && chIdx == 1 
        h5write(tmp4,tmp6,PE,startPEUnAvg,countPEUnAvg);
    end            
elseif strcmpi(DataToLoad,'UnShift')
    if AvgExist
        delete(tmp5); delete(tmp4);
        LF = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
        HF = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
        if LoadPE
            PE = h5read(tmp1,tmp6,startPEUnAvg,countPEUnAvg);
        end      
        LF1 = []; Location = [];
        return      
    elseif UnAvgExist
        delete(tmp4)
        LF = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
        HF = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
        if LoadPE
            PE = h5read(tmp1,tmp6,startPEUnAvg,countPEUnAvg);
        end      
        LF1 = []; Location = [];
        return
    elseif UnShiftExist
        LF = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
        HF = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
        if LoadPE
            PE = h5read(tmp1,tmp6,startPEUnAvg,countPEUnAvg);
        end      
        LF1 = []; Location = [];
        return
    else
        if chIdx == 1
            if ScanPt == 1
                h5create(tmp1,tmp2,sLFUnAvg,'Datatype','double','ChunkSize', csLFUnAvg,'Deflate',5); % Create an h5 file with the LF variable the size of the LF (in double precision)
                h5create(tmp1,tmp3,sHFUnAvg,'Datatype','double','ChunkSize',csHFUnAvg,'Deflate',5);
                % [Fast Time,SlowTime (i.e., Number of Pulses), US Array Elements, Averages, XPt, YPt]
            end
            [HF,LF,PE] = UnAverageReadZeta(Parms,idx,XPt,YPt,LoadPE); % Use the UnAverage reader if the scans have not been averaged           
            tic
            h5write(tmp1,tmp2,LF,startLFUnShift,countLFUnShift); % Write the LF variable to the h5 file
            h5write(tmp1,tmp3,HF,startHFUnShift,countHFUnShift);
            if LoadPE
                h5write(tmp1,tmp6,PE,startPEUnAvg,countPEUnAvg); 
            end
            disp([num2str(toc),' Seconds to Save UnShifted Scan @ X=' num2str(XPt) ' Y=' num2str(YPt)])            
            LF = squeeze(LF(:,:,chIdx));
            HF = squeeze(HF(:,:,:,chIdx));
        else
            LF = h5read(tmp1,tmp2,startLFUnAvg,countLFUnAvg);
            HF = h5read(tmp1,tmp3,startHFUnAvg,countHFUnAvg);
            if LoadPE
                PE = h5read(tmp1,tmp6,startPEUnAvg,countPEUnAvg);
            end                           
        end
    end      
end

if strcmpi(DataToLoad,'UnShift')
    DataLoaded = [Parms.Format.filename2 '_Raw4dUnShifted.h5'];
elseif strcmpi(DataToLoad,'UnAvg')
    DataLoaded = [Parms.Format.filename2 '_Raw4dUnAveraged.h5'];
elseif strcmpi(DataToLoad,'Avg') && isfield(Parms.Format,'Save_Avg') && Parms.Format.Save_Avg
    DataLoaded = [Parms.Format.filename2 '_Raw4d.h5'];
elseif strcmpi(DataToLoad,'Avg') && (~isfield(Parms.Format,'Save_Avg') || ~Parms.Format.Save_Avg)
    DataLoaded = [Parms.Format.filename2 '_Raw4dAveraged.h5'];
elseif strcmpi(DataToLoad,'Filt')
    DataLoaded = [Parms.Format.filename2 '_Filtered4d.h5'];
end
Parms.DataLoaded = DataLoaded;
disp([DataLoaded ' Loaded']);
PE_Exist = exist('PE','var');
if ~PE_Exist
    PE = [];
end
multiWaitbar('CloseAll');
disp([num2str(toc(LoadTimePerPt)),' Seconds to Load all data for chIdx=' num2str(chIdx) ' @ X=' num2str(XPt) ' Y=' num2str(YPt)])                    

end
