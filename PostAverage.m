function [HFAveraged,LFAveraged,LF1] = PostAverage(Parms,LF1,filterCuts,ScanPt,Chan,ManRemoveAvgs,ShiftAverage,AvgToShow,ShowPlots)
%% Arguments
    % Parms (struc)- Parameter file (typically called from bScanParm or equivalent)
    % LF1 (Nx1 double or []) - Low frequency waveform to use for shifting
    % filtercuts (4x1 double) - [Low1 Low2 High1 High2] - Frequency cutoffs for Hamming window
    % ScanPt (double) - Scan Point
    % Chan (double) - Channel Number
    % ManRemoveAvgs (Boolean) - Flag for manually removing averages
    % ShiftAverage (Boolean) - Set to 1
    % AvgToShow (double) - Set to 0
    % Show Pltos (Boolean) - Flag for showing plots during cross-correlation process
%% Outputs
    % HFAveraged (array of double) - Slow Time x Fast Time
    % LFAveraged (array of double) - Amplitude
    % LF1 (Nx1 double) - Low frequency waveform to use for shifting in remainder of scan points
tic
ConversionFactor = (Parms.Daq.HF.PulseRate*Parms.Scan.Duration_ms/1000)/Parms.Daq.LF.Samples; % Determine how many HF samples can be converted to LF Samples
HFChannelNumbers = str2num(strrep(Parms.Daq.HF.Channels,':',' '));
if numel(HFChannelNumbers) > 1
    HFChanns = HFChannelNumbers(1):1:HFChannelNumbers(2); NumHFChanns = numel(HFChanns); % Set the number and names of the HF Channels
else
    HFChanns = HFChannelNumbers; NumHFChanns = 1;
end
LFChannelNumbers = str2num(strrep(Parms.Daq.LF.Channels,':',' '));
if numel(LFChannelNumbers) > 1
    LFChanns = LFChannelNumbers(1):1:LFChannelNumbers(2); NumLFChanns = numel(LFChanns); % Set the number and names of the LF Channels
else
    LFChanns = LFChannelNumbers; NumLFChanns = 1;
end
% Someday include something if HFChans ~= LFChans   
subplotrows = ceil(sqrt(Parms.Scan.Avg));
subplotcolumns = ceil(Parms.Scan.Avg/subplotrows);
%% Read in LF Data
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
blk_size = prod(dsize)*4; % AMAQ: Still don't know why this is multiplied by 4 CP - Has to do with the number of bytes each number is saved as (int32)
fseek(fid,blk_size*(ScanPt-1),'cof'); % Position the marker at the appropriate point for each scan point\
data = fread(fid,fliplr(dsize),'single'); % Flip the orientation of the data matrix
for AvgChanNum = 1:Parms.Scan.Avg*NumLFChanns
    SampsPerAvg = dsize(2)/Parms.Scan.Avg; % Find the number of samples there are in a single average
    SingleAvg(:,AvgChanNum) = data(1:SampsPerAvg); data = data(SampsPerAvg+1:end); % Split the data array into each average
end
fclose(fid);
if ShowPlots
    OriginalPlots = figure(600); %OriginalPlots.Units = 'normalized'; OriginalPlots.OuterPosition = [0 0 1 1]; clf; 
end
for AvgNum = 1:Parms.Scan.Avg
    LFSingle(:,AvgNum) = SingleAvg(:,Chan+(AvgNum-1)*(NumLFChanns)); % Define the LF1 matrix at each scan position as data
    LFSingle(:,AvgNum) = LFSingle(:,AvgNum)-mean(LFSingle(:,AvgNum)); % Remove any DC offset that might be present
    [pks,locs] = findpeaks(LFSingle(:,AvgNum),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(LFSingle(:,AvgNum)),'Annotate','extents');
    [negpks,neglocs] = findpeaks(-LFSingle(:,AvgNum),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(-LFSingle(:,AvgNum)),'Annotate','extents');
    if max(abs(pks)) >= max(abs(negpks))
        peaks = pks; locas = locs;
    else
        peaks = negpks; locas = neglocs;
    end
    AvgExist = exist('AvgsToThrowAway','var');
    if numel(peaks) ~=1          
        if AvgExist ~= 0
            AvgsToThrowAway = cat(2,AvgsToThrowAway,AvgNum);
        else
            AvgsToThrowAway = AvgNum;                
        end
    end
    if ShowPlots  && numel(peaks) == 1
        subplot(subplotrows,subplotcolumns,AvgNum)
        plot(LFSingle(:,AvgNum)); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['ScanPt ' num2str(ScanPt) ' Chann ' num2str(LFChanns(Chan)) ' Avg ' num2str(AvgNum)]); axis tight;  hold on
        plot(locs,pks, ' ro ');  
        plot(neglocs,-negpks,' bo '); hold off
    end
end
    %% Read in HF Data
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
tsize = dsize(1)*dsize(2)*dsize(3)*dsize(4); % AMAQ: I know what dsize(3) represents as the number of burst triggers, but not sure of the others
if tsize ~= length(data) % If tsize is not the same as the data just read in then a missing point exists at this scan point
    data = zeros(tsize,1); % Set all the data to zeros of the appropriate size at Missing Point
    disp(['Missing point ' num2str(ScanPt)]); % Report to the screen which scan point now has zeroes
end
% data = reshape(data,permute(dsize,[4,3,2,1])); % Reshape the data matrix to an 3-D array with fast time on the first dimension and slow time along the last dimension
% HFSingle = permute(data,[2,3,4,2]); % Re-arrange data so that the 2nd and 3rd dimension are switched (i.e., the single unit is moved to the end)

data = reshape(data,[dsize(3),dsize(2),dsize(1),dsize(4)]); % Original [FastTime Channel SlowTime Average] % Reshape the data matrix to an 3-D array with fast time on the first dimension and slow time along the last dimension
% data = reshape(data,permute(dsize,[3 2 1 4])); % Original [FastTime Channel SlowTime Average] % Reshape the data matrix to an 3-D array with fast time on the first dimension and slow time along the last dimension
HFSingle = permute(data,[3,1,4,2]); % Original [SlowTime FastTime Average Channel] % Re-arrange data so that the 2nd and 3rd dimension are switched (i.e., the single unit is moved to the end)
HFSingle = squeeze(HFSingle(:,:,:,Chan));
fclose(fid);            
%% Throw out averages without a peak in them 
% In the Future include a section here that will do wavelet-based peak detection to rule in signals
% For now, just ask which averages the user should keep
if ManRemoveAvgs
    USAvgsToThrowAway = inputdlg('Select Which Averages to Throw Away (Separate by Space)','Averages to Throw');
    AvgsToThrowAway = cat(2, AvgsToThrowAway,str2num(USAvgsToThrowAway{:})); % Convert to a double
end
AvgExist = exist('AvgsToThrowAway','var');
if AvgExist ~= 0
    LFSingle(:,AvgsToThrowAway) = []; LFSingle = squeeze(LFSingle);  % Set Averages that the user does not want to include in averaging to null and recreate size of arrzy
    HFSingle(:,:,AvgsToThrowAway) = []; HFSingle = squeeze(HFSingle); % Repeat for the HF
end
%% Perform cross-correlation
NumAvgs = size(LFSingle,2); % Determine the new number of averages after others have been thrown away
subplotrows2 = ceil(sqrt(NumAvgs)); % Determine an appropriate subplot array for plotting
subplotcolumns2 = ceil((NumAvgs)/subplotrows2);     
if ShowPlots
    gcf; clf; % Using the figure on which the original figures were plotted
end
for AvgNum = 1:NumAvgs
    if ShowPlots
        subplot(subplotrows2,subplotcolumns2,AvgNum); % Plot and add labels and titles to each plot
        plot(LFSingle(:,AvgNum)); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['ScanPt ' num2str(ScanPt) ' Chann ' num2str(LFChanns(Chan)) ' Avg ' num2str(AvgNum)]); axis tight; 
    end
end
if ShiftAverage
    if isempty(LF1)
        LF1 = LFSingle(:,1); % Set the first LF1 signal equal to the first average
    end
    if ShowPlots
        figure(700); clf; figure(800); clf; % Open two new figures for plotting cross-correlations and shifted signals
    end
    LFAveraged = zeros(size(LF1)); % Sum the LF waveform
    if size(HFSingle,3) >= 1
        HFAveraged = zeros(size(HFSingle(:,:,1))); % Sum the HF waveform 
    else
        HFAveraged = zeros(size(HFSingle));
    end
    padLF = 40;
    LF1 = LF1';
    LF1 = padarray(LF1',[padLF,0],'replicate','both')';
    LF1 = 2*real(fft_filt2D_new(LF1',[],Parms.Daq.LF.Rate_hz*1000,filterCuts,0))';
    LF1 = LF1(:,1+padLF:end-padLF)';
    for AvgNum = 1:NumAvgs
        LFShifted = zeros(size(LF1)); % Create an array of zeros to put the shifted LF into     
        LF2 = LFSingle(:,AvgNum); % Set the LF2 signal equal to each average (including the first)
        LFFilt = LF2';
        LFFilt = padarray(LFFilt',[padLF,0],'replicate','both')';
        LFFilt = 2*real(fft_filt2D_new(LFFilt',[],Parms.Daq.LF.Rate_hz*1000,filterCuts,0))';
        LFFilt = LFFilt(:,1+padLF:end-padLF)';
        [cors, lags] = xcorr(abs(LF1),abs(LFFilt)); % Calculate the lags and correlation values between the LF1 and LF2 signal
        [~,I] = max(abs(cors)); % Find the index of the maximum absolute correlation
    %     if strcmp(Parms.Stim.Waveform,'Arb') || strcmp(Parms.Stim.Waveform,'ECG') || strcmp(Parms.Stim.Waveform,'Sin') % Take out Sin when you actually start saving this correctly 
    %         [corpks,corlocs] = findpeaks(cors,'WidthReference','halfprom','MaxPeakWidth',0.1*numel(cors),'MinPeakProminence',0.1*max(abs(cors)),'Threshold',5e-7,'Annotate','extents');
    %         if numel(corpks) ~=1
    %             lagDiff = lags(I); % Find the lag that corresponds to that maximum
    %         else
    %             lagDiff = lags(corlocs);
    %         end
    %     else
            lagDiff = lags(I);
    %     end
        ConvFactor = ceil(lagDiff*ConversionFactor); % Calculate how many HF samples this lag corresponds to
        % Plot the correlations
        if ShowPlots
            figure(700);
            subplot(subplotrows2,subplotcolumns2,AvgNum)
            plot(lags,cors); xlabel('Lags'); ylabel('Correlation Value'); title(['Correlation between 1&' num2str(AvgNum)]); axis tight;      
            hold on
    %         plot(corlocs+min(lags),corpks, ' ro ')
            hold off
        end
        %% Shift Low and High Frequency Peaks
        if lagDiff > 0 % If the peak in LF2 occurs prior to the peak in LF1, then shift to the right
            lagDiff = lagDiff;
            LFShifted(lagDiff:end) = LF2(1:size(LFShifted(lagDiff:end))); 
        elseif lagDiff < 0 % If the peak in LF2 occurs after the peak in LF1, then shift to the left
            lagDiff = -lagDiff;
            LFShifted(1:end-lagDiff) = LF2(lagDiff+1:end); 
        else 
    %                     lagDiff = 0; % If the peak in LF2 occurs at the same point as the peak in LF1, then use the whole waveform
            lagDiff = 1;
            LFShifted(lagDiff:end) = LF2(1:size(LFShifted(lagDiff:end))); 
        end
    %                 LFShifted = circshift(LF2,lagDiff);             
        HF2 = HFSingle(:,:,AvgNum); % Create a dummy variable HF2 for shifting each average of data
        HFShifted = zeros(size(HF2)); % Create an array of zeros for shifted HF peaks (same as LF)
        % Repeat all shifting as in LF
        if ConvFactor > 0
            ConvFactor = ConvFactor;
            HFShifted(ConvFactor:end,:) = HF2(1:size(HFShifted(ConvFactor:end,:),1),:);
        elseif ConvFactor < 0
            ConvFactor = -ConvFactor;
            HFShifted(1:end-ConvFactor,:) = HF2(ConvFactor+1:end,:);
        else 
    %                     ConvFactor = 0;
            ConvFactor = 1;
            HFShifted(ConvFactor:end,:) = HF2(1:size(HFShifted(ConvFactor:end,:),1),:);
        end
    %                 HFShifted = circshift(HF2,lagDiff,1);                             
        % Plot newly shifted LF waveforms
        if ShowPlots
            figure(800);
            subplot(subplotrows2,subplotcolumns2,AvgNum)
            plot(LFShifted); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['ScanPt ' num2str(ScanPt) ' Chann ' num2str(LFChanns(Chan)) ' Avg ' num2str(AvgNum)]); axis tight;   
        end
        LFAveraged = LFAveraged+LFShifted; % Sum the LF waveform
        HFAveraged = HFAveraged+HFShifted; % Sum the HF waveform

    end
    LFAveraged = LFAveraged./NumAvgs; % Average the LF waveform
    HFAveraged = HFAveraged./NumAvgs; % Average the HF waveform
    % Plot the averaged LF Waveform
    if ShowPlots
        figure(900+Chan+ScanPt); clf;
        plot(LFAveraged); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['Averaged ScanPt ' num2str(ScanPt) ' LF Chann ' num2str(LFChanns(Chan))]); axis tight; % Eventually will be LFAveraged(:,XPt,YPt,Chan)     

    end
else
    LFAveraged = LFSingle(:,AvgToShow);
    HFAveraged = HFSingle(:,:,AvgToShow);
    LF1 = [];
end
HFAveraged = permute(HFAveraged,[2 1]);
disp([num2str(toc),' Seconds to Average Scan'])
end
