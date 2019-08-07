function [HFSingle, LFSingle, PESingle] = UnAverageReadZeta(Parms, ScanPt, XPt,YPt, LoadPE)
%% UnAverageRead.m - Reader for data from NI 5105 and 6339 that does not average while collecting
% Written by Alex Alvarez
% Last Updated 04/19/19
% Arguments
%   Parms (struc)       - Structural array containing all scan parameters
%   ScanPt (double)     - Index of the X/Y Scan
%   XPt (double)        - Index of the X Scan
%   YPt (double)        - Index of the Y Scan
%   LoadPE (Boolean)    - Flag to LoadPE or not
% Output
%   LFSingle (3D Array) - New (SlowTime x Average x Channel)
%   HFSingle (4D Array) - New (FastTime x SlowTime x Average x Channel)
%   PESingle (4D Array of Doubles) -  (FastTime x Slow Time (i.e., Pulse Number) x Element x Average]
%% TOC
% @00 - High-Level Flags
% @01 - Definitions
% @02 - Read LF Data
% @03 - Read HF Data
% @04 - Read PE Data
%% @00 - High-Level Flags
%% @01 - Initial Definitions
% If there is more than one channel, then define LFChanns and NumLFChanns
% accordingly
ReadTime = tic;
LFChannelNumbers = str2num(strrep(Parms.Daq.LF.Channels,':',' '));
if numel(LFChannelNumbers) > 1
    LFChanns = LFChannelNumbers(1):1:LFChannelNumbers(2); NumLFChanns = numel(LFChanns); % Set the number and names of the LF Channels
else
    LFChanns = LFChannelNumbers; NumLFChanns = 1;
end
%% @02 - Read in LF Data
lffile = strcat(Parms.Format.filename2,'_LF.dat'); % Define the low frequency average file
fid = fopen(lffile); % Open the file ID for the file
if fid < 0 % If an fid cannot be created
    error('Cannot open LF data file') % Report an error message to the screen
end
n = fread(fid,1,'int32'); % Read the number of characters to include in dsize
sztype = fgets(fid,n); % Read in new line with the number of characters defined by n
if(strcmpi(sztype(1:5),'ucsdi') == 0) % If this new line does not contain ucsdi
    fclose(fid); % Then close the file
    return; % And return to the workspace
end
Vers = fread(fid,1,'int32'); % Move the pointer past the date version value
n = fread(fid,1,'int32'); % Read the number of characters to include in dsize
dsize = fread(fid,[1,n],'int32'); % Define the size of the LF array based on the number of characters read
blk_size = prod(dsize)*4; % AMAQ: Still don't know why this is multiplied by 4
fseek(fid,blk_size*(ScanPt-1),'cof'); % Position the marker at the appropriate point for each scan point\
% % data = fread(fid,fliplr(dsize),'single'); % Flip the orientation of the data matrix
data = fread(fid,fliplr(dsize),'single'); % Flip the orientation of the data matrix
if isempty(data) || size(data,1) ~= dsize(2) || size(data,2) ~= dsize(1)
    data = zeros(fliplr(dsize));
end
data = reshape(data,[prod(dsize),1]);
for AvgChanNum = 1:Parms.Scan.Avg*NumLFChanns
    SampsPerAvg = dsize(2)/Parms.Scan.Avg; % Find the number of samples there are in a single average
    SingleAvg(:,AvgChanNum) = data(1:SampsPerAvg); data = data(SampsPerAvg+1:end); % Split the data array into each average
end
for Chan = 1:NumLFChanns
    for AvgNum = 1:Parms.Scan.Avg % Define the
        LFSingle(:,AvgNum,Chan) = SingleAvg(:,Chan+(AvgNum-1)*(NumLFChanns)); % Define the LF1 matrix at each scan position as data
    end
end

fclose(fid);
%% @03 - Read in HF Data
hffile = strcat(Parms.Format.filename2,'_HF.dat'); % Define the name of the HF file
if exist(hffile,'file') % If this file exists
    fid = fopen(hffile,'rb'); % Define the file id of the file and allow it to be readable as a non-text file
end
if fid < 0 % If the file does not exist
    error('Cannot open HF data file'); % Output error message
end
n = fread(fid,1,'int32'); % Define the number of characters as the first integer value of the file
sztype = fgets(fid,n); % Read in new line with the number of characters defined by n
if(strcmpi(sztype(1:5),'ucsdi') == 0) % If this new line does not contain ucsdi
    fclose(fid); % Then close the file
    return; % And return to the workspace
end
Vers = fread(fid,1,'int32'); % Move the pointer past the date version value
n = fread(fid,1,'int32'); % Determine the 3D size to read the array into 
dsize = fread(fid,[1,n],'int32'); % [SlowTime Channel FastTime Average] Read in the dimensions of the array (Number of time points x Number of dimensions x Number of samples)
blk_size = prod(dsize)*4; % AMAQ: Clearly the product of the dimensions by why multiplied by 4 What is the blk_size supposed to represent
fseek(fid,blk_size*(ScanPt-1),'cof'); % Move the file marker in fid to blk_size * the ScanPt-1 starting from the current position in the file so as to read one block of data for each scan
data = fread(fid,[blk_size/4,1],'single'); % Read in a vector reprsenting the total number of fast time samples * slow time points into data array with single precision
tsize = dsize(1)*dsize(2)*dsize(3)*dsize(4); % Frames x NumofHFChans x HFSamples x NumofAvg
% if tsize ~= length(data) % If tsize is not the same as the data just read in then a missing point exists at this scan point
%     data = zeros(tsize,1); % Set all the data to zeros of the appropriate size at Missing Point
%     disp(['Missing point ' num2str(ScanPt)]); % Report to the screen which scan point now has zeroes
% end
if isempty(data) || size(data,1) ~= tsize 
    data2 = zeros(dsize(3),dsize(2),dsize(1)*dsize(4));
else
    data2 = reshape(data,[dsize(3),dsize(2),dsize(1)*dsize(4)]);
end
    % data2 = reshape(data,[dsize(3),dsize(2),dsize(1)*dsize(4)]);

for AvgChanNum = 1:Parms.Scan.Avg
    SampsPerAvg = dsize(1); % Find the number of samples there are in a single average
    SingleHFAvg(:,:,:,AvgChanNum) = data2(:,:,1:SampsPerAvg); data2 = data2(:,:,SampsPerAvg+1:end); % Split the data array into each average
end
HFSingle = permute(SingleHFAvg,[1,3,4,2]); % Original [FastTime SlowTime Average Channel] % Re-arrange data so that the 2nd and 3rd dimension are switched (i.e., the single unit is moved to the end)
fclose(fid); % Close the fid
%% @04 - Load PE Data
if LoadPE 
    cd ..\PEData % change to the appropriate path location for PE data
    DisStruc = load([Parms.Format.filename2 '_PEParm.mat']); % Load the path
    bScanParm = DisStruc.bScanParm; % [struc] Defines the parameters used for the scans
    Trans = DisStruc.Trans; % [struc] Defines the transducer structure used in scan
    Receive = DisStruc.Rcv; % [struc] Defines the receive structure used in scan
%     Resource = DisStruc.Resource; % [struc] Defines the resource structure used in scan
    ActFrames = bScanParm.ActualSammies*bScanParm.NumyFrames;
    
    if ActFrames <= 4096*128
        Frames = 2;
    else
        [~,Frames] = min(ActFrames-(4096*128*linspace(1,20,20)));
    end    
    % size(RF_Data) = [Samples x Channels x Frames x Averages x XPt x YPt]
    AcqPerFrame = bScanParm.NumyFrames; % [double] Number of acquisition pulses per frame
    Acqs = Receive(AcqPerFrame).endSample-Receive(1).startSample+1;
%     Frames = Resource.RcvBuffer(1).numFrames;
    NumChanns = bScanParm.MModeCols;
    countPE =  [Acqs,NumChanns,Frames,Parms.Scan.Avg,1,1];
    startPE =  [1,1,1,1,XPt,YPt];
    RF_Data = h5read([bScanParm.Format.filename2 '_MModePerPulse.h5'],'/PEMMode',startPE,countPE); % [doubles] array of PE data saved in h5
    AcqSize = Receive(1).endSample - Receive(1).startSample + 1; % [double] Size of each pulse acquisition
    DataSize = size(RF_Data); % [doubles] Vector of sizes for PEMMode 
    %  Reshape A into the appropriate dimensions based on size of A; E
    %  is an array of doubles
    if numel(DataSize) == 4
        ReshapeRF = reshape(RF_Data,[AcqSize, AcqPerFrame, DataSize(2), DataSize(3), DataSize(4)]); % [Samples per acq, PulseNum, Channels, Frames, Averages]
    elseif numel(DataSize) == 5 && bScanParm.Scan.Xpt > 1
        ReshapeRF = reshape(RF_Data,[AcqSize, AcqPerFrame, DataSize(2), DataSize(3), DataSize(4), DataSize(5)]);% [Samples per acq, PulseNum, Channels, Frames, Averages, XPts]
        ReshapeRF = squeeze(ReshapeRF(:,:,:,:,:,XPt));
    else
        ReshapeRF = reshape(RF_Data,[AcqSize, AcqPerFrame, DataSize(2), DataSize(3), DataSize(4), DataSize(5), DataSize(6)]);% [Samples per acq, PulseNum, Channels, Frames, Averages, XPts, YPts]
        % size(ReshapeRF) = [Fast Time x Slow Time x Channels x Frames x Averages]
        ReshapeRF = squeeze(ReshapeRF(:,:,:,:,:,XPt,YPt));
    end
    % If there is a single frame acquired then squeeze E; otherwise
    % create an array out of each of the frames
    if ActFrames <= 4096*128
        ReshapeRF = squeeze(ReshapeRF(:,:,:,1,:));
    else
        for ReshapeSize = 1:size(ReshapeRF,4)
            if ReshapeSize == 1
                ReshapeRF = squeeze(ReshapeRF(:,:,:,ReshapeSize,:)); % [Samples per acq, PulseNum, Channels, Averages]
            else
                ReshapeRF = cat(2,squeeze(ReshapeRF(:,:,:,ReshapeSize,:)),ReshapeRF); % [Samples per acq, PulseNum, Channels, Averages]
            end
        end
    end
    % size(ReshapeRF) = [Fast Time x Slow Time x Channels x Averages]
    NumChanns = size(ReshapeRF,3); % [double] Number of channels that have RF data

    if bScanParm.twoAcq
        % If there were an even number of pulses taken then use all the
        % pulses; otherwise use all the pulses except the last one (this
        % should only apply to MultiAcq)
        if mod(size(ReshapeRF,2),2) == 0 
            SqueezeSize = size(ReshapeRF,2); % [double] The size of the squeezed array
        else
            SqueezeSize = size(ReshapeRF,2)-1;
        end           
        % size(SynthApRF) = [Samples per acq, PulseNum/2, Channels, Averages]
        SynthApRF = zeros(size(ReshapeRF,1),SqueezeSize/2,size(ReshapeRF,3),size(ReshapeRF,4)); % [doubles] An array that will be filled when combining two consecutive pulses with half the aperture
        ActPulseNum = 0; SynthPulseNum = 0;
        % Loop over h and j to create a new array I that will combine
        % the two pulses of a synthetic aperture into a single pulse
        while ActPulseNum < SqueezeSize          
            ActPulseNum = ActPulseNum+1;
            if mod(ActPulseNum,2) ~= 0
                SynthPulseNum = SynthPulseNum+1;
                SynthApRF(:,SynthPulseNum,1:NumChanns/2,:) = ReshapeRF(:,ActPulseNum,1:NumChanns/2,:);
            else
                SynthApRF(:,SynthPulseNum,NumChanns/2+1:NumChanns,:) = ReshapeRF(:,ActPulseNum,NumChanns/2+1:end,:);
            end
        end
    else % If not doing two acquisition mode with synthetic aperture, then just define I as the full set of data          
        SynthApRF = ReshapeRF; 
    end
    % size(PESingle) = [Fast Time, Slow Time, Num Elements, Averages]
    PESingle = zeros(size(SynthApRF,1),size(SynthApRF,2),size(Trans.ConnectorES,1),size(SynthApRF,4)); % [doubles] Create an array to fill with data organized by element and not by channel
    % Iterate over j to create array L that will contain the
    % element-mapped data whereas I was the channel-mapped data
    for SynthPulseNum = 1:NumChanns
       k = find(SynthPulseNum == Trans.ConnectorES);
       if ~isempty(k)
            PESingle(:,:,k,:) = SynthApRF(:,:,SynthPulseNum,:);
       end 
    end             
    if ~bScanParm.twoAcq
        if ~bScanParm.twoD
            FindCenters = (Trans.numelements-64)/2;
            PESingle(:,:,end-FindCenters+1:end,:) = []; PESingle(:,:,1:FindCenters,:) = [];
            PESingle = squeeze(PESingle);
        else
            if strcmp(Trans.name(1:4),'H235')
                % Set to the center 60 elements in the array
                PESingle(:,:,106:end,:) = []; PESingle(:,:,88:93,:) = []; PESingle(:,:,70:75,:) = [];  
                PESingle(:,:,52:57,:) = []; PESingle(:,:,34:39,:) = []; PESingle(:,:,1:21,:) = [];
                PESingle = squeeze(PESingle);
            elseif strcmp(Trans.name(1:4),'H247')
                PESingle(:,:,105:end,:) = [];
                PESingle(:,:,98:99,:) = [];
                PESingle(:,:,91:92,:) = [];
                PESingle(:,:,84:85,:) = [];
                PESingle(:,:,77:78,:) = [];
                PESingle(:,:,70:71,:) = [];
                PESingle(:,:,63:64,:) = [];
                PESingle(:,:,56:57,:) = [];
                PESingle(:,:,49:50,:) = [];
                PESingle(:,:,42:43,:) = [];
                PESingle(:,:,35:36,:) = [];
                PESingle(:,:,28:29,:) = [];
                PESingle(:,:,1:22,:) = [];
                PESingle = squeeze(PESingle);
            end
        end
    end    
    cd ..\ExpData    
else
    PESingle = [];
end
disp([num2str(toc(ReadTime)),' Seconds to Read Scan @ X=' num2str(XPt) ' Y=' num2str(YPt)])
