function [HFAveraged,LFAveraged,LF1,NumThrown,NumSplit] = PostAverageDelta(Parms,LF1,SingleAvg,HFSingle,filterCuts,ScanPt,Chan,Manually,ShowPlots)
%% Changes to be Included in Beta
% 1. Reader will be separated from shifter/averager (Done through
% UnAverageRead)
% 2. Template will be shifted to center (integrated in PostAverageBeta)
% 3. Averages with two peaks will be split into two averages (integrated in PostAverageBeta)
%% Arguments
    % Parms (struc)- Parameter file (typically called from bScanParm or equivalent)
    % LF1 (Nx1 double or []) - Low frequency waveform to use for shifting
    % SingleAvg (Nx1 double) - Low frequency data set from all averages and one channel at one scan point (defined by Chan and ScanPt)
    % HFSingle (NxMxP double) - High frequency data set from all averages and one channel at one scan point (defined by Chan and ScanPt)   
    % filtercuts (4x1 double) - [Low1 Low2 High1 High2] - Frequency cutoffs for Hanning window
    % ScanPt (double) - Scan Point
    % Chan (double) - Channel Number
    % Manually (Boolean) - Flag for manually removing averages
    % Show Plots (Boolean) - Flag for showing plots during cross-correlation process
%% Outputs
    % HFAveraged (array of double) - Slow Time x Fast Time
    % LFAveraged (array of double) - Amplitude
    % HFShiftedMat (array of double) - Includes all averages of HF data before averaging together
    % LFShiftedMat (array of double) - Same for LF
    % LF1 (Nx1 double) - Low frequency waveform to use for shifting in remainder of scan points
    % NumThrown (double) - Number of Averages thrown during the averaging

%% Initial definitions
tic % Start a clock to see how long it takes
ConversionFactor = (Parms.Daq.HF.PulseRate*Parms.Scan.Duration_ms/1000)/Parms.Daq.LF.Samples; % Determine how many HF samples can be converted to LF Samples
HFChannelNumbers = str2num(strrep(Parms.Daq.HF.Channels,':',' ')); % Determine the HF Channel Numbers
if numel(HFChannelNumbers) > 1
    HFChanns = HFChannelNumbers(1):1:HFChannelNumbers(2); NumHFChanns = numel(HFChanns); % Set the number and names of the HF Channels
else
    HFChanns = HFChannelNumbers; NumHFChanns = 1;
end
LFChannelNumbers = str2num(strrep(Parms.Daq.LF.Channels,':',' ')); % Determine the LF channel numbers used
if numel(LFChannelNumbers) > 1
    LFChanns = LFChannelNumbers(1):1:LFChannelNumbers(2); NumLFChanns = numel(LFChanns); % Set the number and names of the LF Channels
else
    LFChanns = LFChannelNumbers; NumLFChanns = 1;
end
% Someday include something if HFChans ~= LFChans   
subplotrows = ceil(sqrt(Parms.Scan.Avg)); % Determine the number of rows needed for an appropriate subplot array for the number of averages
subplotcolumns = ceil(Parms.Scan.Avg/subplotrows); % Do the same with the columns
%% Read in LF Data
if ShowPlots % Create a figure if the flag for ShowPlots is set to 1
    OriginalPlots = figure(600); OriginalPlots.Units = 'normalized'; OriginalPlots.OuterPosition = [0 0 1 1]; clf; % Set the figure to take the entirety of the screen
end
for AvgNum = 2:Parms.Scan.Avg % Define the 
    LFSingle(:,AvgNum) = SingleAvg(:,AvgNum)-mean(SingleAvg(:,AvgNum)); % Remove any DC offset that might be present
    [pks,locs] = findpeaks(LFSingle(:,AvgNum),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(LFSingle(:,AvgNum)),'Annotate','extents');
    [negpks,neglocs] = findpeaks(-LFSingle(:,AvgNum),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(-LFSingle(:,AvgNum)),'Annotate','extents');
    if max(abs(pks)) >= max(abs(negpks))
        peaks = pks; locas = locs;
    else
        peaks = negpks; locas = neglocs;
    end
    if ~Manually
        AvgExist = exist('AvgsToThrowAway','var');
        SplitExist = exist('AvgsToSplit','var');
        if numel(peaks) == 0          
            if AvgExist ~= 0
                AvgsToThrowAway = cat(2,AvgsToThrowAway,AvgNum);
            else
                AvgsToThrowAway = AvgNum;                
            end
        elseif numel(peaks) == 2
            if SplitExist ~= 0
                AvgsToSplit = cat(2,AvgsToSplit,AvgNum);
            else
                AvgsToSplit = AvgNum;                
            end        
        end
    end
    if ShowPlots  % && numel(peaks) == 1
        subplot(subplotrows,subplotcolumns,AvgNum)
        plot(LFSingle(:,AvgNum)); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['Pt' num2str(ScanPt) ' Chan' num2str(LFChanns(Chan)) ' Avg' num2str(AvgNum)]); axis tight;  hold on
        plot(locs,pks, ' ro ');  
        plot(neglocs,-negpks,' bo '); hold off
    end
end
%% Read in HF Data
%% Throw out averages without a peak in them 
% In the Future include a section here that will do wavelet-based peak detection to rule in signals
% For now, just ask which averages the user should keep
if Manually
    USAvgsToThrowAway = inputdlg('Select Which Averages to Throw Away (Separate by Space)','Averages to Throw');
    AvgsToThrowAway = str2num(USAvgsToThrowAway{:}); %cat(2, AvgsToThrowAway,str2num(USAvgsToThrowAway{:})); % Convert to a double
    USAvgsToSplit = inputdlg('Select Which Averages to Split in Two (Separate by Space)','Averages to Split');
    AvgsToSplit = str2num(USAvgsToSplit{:}); 
end
AvgExist = exist('AvgsToThrowAway','var');
SplitExist = exist('AvgsToSplit','var');
if SplitExist ~= 0
    NumSplit = numel(AvgsToSplit);
    for q = 1:NumSplit
        SplitAvg = AvgsToSplit(q);
        [pks,locs] = findpeaks(LFSingle(:,SplitAvg),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(LFSingle(:,SplitAvg)),'Annotate','extents');
        [negpks,neglocs] = findpeaks(-LFSingle(:,SplitAvg),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(-LFSingle(:,SplitAvg)),'Annotate','extents');        
        if (max(abs(pks)) >= max(abs(negpks)))
            peaks = pks; locas = locs;
        elseif isempty(negpks)
            peaks = pks; locas = locs;
        else
            peaks = negpks; locas = neglocs;
        end
        
        SplitPoint = round((locas(2)-locas(1))/2); SplitPointHF = round(SplitPoint * ConversionFactor);
        FirstHalf = LFSingle(1:SplitPoint,SplitAvg); FirstHalfHF = HFSingle(1:SplitPointHF,:,SplitAvg);
        SecondHalf = LFSingle(SplitPoint+1:end,SplitAvg); SecondHalfHF = HFSingle(SplitPointHF+1:end,:,SplitAvg);
        PadFirstHalf = size(LFSingle(:,SplitAvg),1)-size(FirstHalf,1);
        LFSingle(:,SplitAvg) = padarray(FirstHalf,[PadFirstHalf,0],'replicate','post');
        PadSecondHalf = size(LFSingle(:,SplitAvg),1)-size(SecondHalf,1);
        LFSingle(:,end+1) = padarray(SecondHalf,[PadSecondHalf,0],'replicate','pre');
        PadFirstHalfHF = size(HFSingle(:,:,SplitAvg),1)-size(FirstHalfHF,1);
        HFSingle(:,:,SplitAvg) = padarray(FirstHalfHF,[PadFirstHalfHF,0,0],'replicate','post');
        PadSecondHalfHF = size(HFSingle(:,:,SplitAvg),1)-size(SecondHalfHF,1);
        HFSingle(:,:,end+1) = padarray(SecondHalfHF,[PadSecondHalfHF,0,0],'replicate','pre');
    end
end
if AvgExist ~= 0
    NumThrown = numel(AvgsToThrowAway);
    LFSingle(:,AvgsToThrowAway) = []; LFSingle = squeeze(LFSingle);  % Set Averages that the user does not want to include in averaging to null and recreate size of arrzy
    HFSingle(:,:,AvgsToThrowAway) = []; HFSingle = squeeze(HFSingle); % Repeat for the HF
else
    NumThrown = 0;
end
%% Perform cross-correlation
NumAvgs = size(LFSingle,2); % Determine the new number of averages after others have been thrown away
subplotrows2 = ceil(sqrt(NumAvgs)); % Determine an appropriate subplot array for plotting
subplotcolumns2 = ceil((NumAvgs)/subplotrows2);     

if ShowPlots
    gcf; clf; % Using the figure on which the original figures were plotted
    pause(0.01) % Wait for 1/100 of a second
end
for AvgNum = 1:NumAvgs % Iterate over all averages to plot each LF average
    if ShowPlots
        subplot(subplotrows2,subplotcolumns2,AvgNum); % Plot and add labels and titles to each plot
        plot(LFSingle(:,AvgNum)); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['Pt' num2str(ScanPt) ' Chan' num2str(LFChanns(Chan)) ' Avg' num2str(AvgNum)]); axis tight; 
    end
end
if isempty(LF1)
    LF1 = LFSingle(:,1); % Set the first LF1 signal equal to the first average
    [pks,locs] = findpeaks(LF1,'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(LF1),'Annotate','extents');
    [negpks,neglocs] = findpeaks(-LF1,'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(-LF1),'Annotate','extents');

    if (max(abs(pks)) >= max(abs(negpks)))
        peaks = pks; locas = locs;
    elseif isempty(negpks)
        peaks = pks; locas = locs;
    else
        peaks = negpks; locas = neglocs;
    end
    [~, maxloc] = max(peaks(:));
    location = locas(maxloc);
    CenterLag = location-round(Parms.Daq.LF.Samples/2);
    if CenterLag < 0
        CenterLag = -CenterLag;
        LF1Shifted = LF1(1:end-CenterLag,:);
        PadLF1 = size(LF1,1)-size(LF1Shifted,1);
        LF1Shifted = padarray(LF1Shifted,[PadLF1,0],'replicate','pre');            
    elseif CenterLag > 0
        CenterLag = CenterLag;
        LF1Shifted = LF1(CenterLag+1:end,:);
        PadLF1 = size(LF1,1)-size(LF1Shifted,1);
        LF1Shifted = padarray(LF1Shifted,[PadLF1,0],'replicate','post');
    else 
        CenterLag = 1;
        LF1Shifted = LF1(1:end-CenterLag,:);
        PadLF1 = size(LF1,1)-size(LF1Shifted,1);
        LF1Shifted = padarray(LF1Shifted,[PadLF1,0],'replicate','pre');
    end
    LF1 = LF1Shifted;
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
for AvgNum = 2:NumAvgs
    LFShifted = zeros(size(LF1)); % Create an array of zeros to put the shifted LF into     
    LF2 = LFSingle(:,AvgNum); % Set the LF2 signal equal to each average (including the first)
    LFFilt = LF2';
    LFFilt = padarray(LFFilt',[padLF,0],'replicate','both')';
    LFFilt = 2*real(fft_filt2D_new(LFFilt',[],Parms.Daq.LF.Rate_hz*1000,filterCuts,0))';
    LFFilt = LFFilt(:,1+padLF:end-padLF)';
    [cors, lags] = xcorr(abs(LF1),abs(LFFilt)); % Calculate the lags and correlation values between the LF1 and LF2 signal
    [~,I] = max(abs(cors)); % Find the index of the maximum absolute correlation
    lagDiff = lags(I);
    ConvFactor = round(lagDiff*ConversionFactor); % Calculate how many HF samples this lag corresponds to
    % Plot the correlations
    if ShowPlots
        figure(700);
        subplot(subplotrows2,subplotcolumns2,AvgNum)
        plot(lags,cors); xlabel('Lags'); ylabel('Correlation Value'); title(['Correlation between 1&' num2str(AvgNum)]); axis tight;      
        hold on
    end
    %% Shift Low and High Frequency Peaks
    if lagDiff > 0 % If the peak in LF2 occurs prior to the peak in LF1, then shift to the right
        lagDiff = lagDiff;
%             LFShifted(lagDiff:end) = LF2(1:size(LFShifted(lagDiff:end))); 
        LFShifted = LF2(1:end-lagDiff);
        PadLF = numel(LF2)-numel(LFShifted);
        LFShifted = padarray(LFShifted,PadLF,'replicate','pre');
    elseif lagDiff < 0 % If the peak in LF2 occurs after the peak in LF1, then shift to the left
        lagDiff = -lagDiff;
%             LFShifted(1:end-lagDiff) = LF2(lagDiff+1:end); 
        LFShifted = LF2(lagDiff+1:end);
        PadLF = numel(LF2)-numel(LFShifted);
        LFShifted = padarray(LFShifted,PadLF,'replicate','post');            

    else 
%                         lagDiff = 0; % If the peak in LF2 occurs at the same point as the peak in LF1, then use the whole waveform
        lagDiff = 1;
%             LFShifted(lagDiff:end) = LF2(1:size(LFShifted(lagDiff:end)));
        LFShifted = LF2;
    %    PadLF = numel(LF2)-numel(LFShifted);
     %   LFShifted = padarray(LFShifted,PadLF,'replicate','pre');
    end
%         LFShifted(1:end-lagDiff) = LF2(lagDiff+1:end);        
%                         LFShifted = circshift(LF2,lagDiff,1);
    HF2 = HFSingle(:,:,AvgNum); % Create a dummy variable HF2 for shifting each average of data
%         HFShifted = zeros(size(HF2)); % Create an array of zeros for shifted HF peaks (same as LF)
    % Repeat all shifting as in LF
    if ConvFactor > 0
        ConvFactor = ConvFactor;
%             HFShifted(ConvFactor:end,:) = HF2(1:size(HFShifted(ConvFactor:end,:),1),:);
        HFShifted = HF2(1:end-ConvFactor,:);
        PadHF = size(HF2,1)-size(HFShifted,1);
        HFShifted = padarray(HFShifted,[PadHF,0],'replicate','pre');            
    elseif ConvFactor < 0
        ConvFactor = -ConvFactor;
%             HFShifted(1:end-ConvFactor,:) = HF2(ConvFactor+1:end,:); %             zeroPadding
        HFShifted = HF2(ConvFactor+1:end,:);
        PadHF = size(HF2,1)-size(HFShifted,1);
        HFShifted = padarray(HFShifted,[PadHF,0],'replicate','post');
    else 
%                         ConvFactor = 0;
        ConvFactor = 1;
%             HFShifted(ConvFactor:end,:) = HF2(1:size(HFShifted(ConvFactor:end,:),1),:);
        HFShifted = HF2;
       % PadHF = size(HF2,1)-size(HFShifted,1);
      %  HFShifted = padarray(HFShifted,[PadHF,0],'replicate','pre');
    end
%         HFShifted(1:end-ConvFactor,:) = HF2(ConvFactor+1:end,:);
%                     HFShifted = circshift(HF2,ConvFactor,2);                             
    % Plot newly shifted LF waveforms
    if ShowPlots
        figure(800);
        subplot(subplotrows2,subplotcolumns2,AvgNum)
        plot(LFShifted); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['Pt' num2str(ScanPt) ' Chan' num2str(LFChanns(Chan)) ' Avg' num2str(AvgNum)]); axis tight;   
    end
    LFAveraged = LFAveraged+LFShifted; % Sum the LF waveform
    HFAveraged = HFAveraged+HFShifted; % Sum the HF waveform
end
LFAveraged = LFAveraged./NumAvgs; % Average the LF waveform
HFAveraged = HFAveraged./NumAvgs; % Average the HF waveform

% Plot the averaged LF Waveform
figure(900+Chan);
plot(LFAveraged); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['Averaged  ' ' LF Chann ' num2str(LFChanns(Chan))]); axis tight; % Eventually will be LFAveraged(:,XPt,YPt,Chan)     
hold on
HFAveraged = permute(HFAveraged,[2 1]);
disp([num2str(toc),' Seconds to Average Scan'])
end
function [filtSig2,datafft,datafft2,mask2D,freqArray]=fft_filt2D_new(data,fRng,freq,filterCuts,plotit)
%% fft_filt2D_new will take a 1D or 2D matrix and filter matrix--make sure
%%                time is represented by first dimension.
%% in:
%% data = real data in
%% fRng = subset data in points (indices) between fRng(1) and fRng(2), 0 = all
%% freq = sample rate (same dimensions as filter Cuts)
%% filterCuts = 4 numbers in bracket denoting edges of hamming filter (0
%%              outside of filter), e.g., [2 5 8 11] = hamming rise fom 2
%%              to 5, flat top, hanning fall 8 to 11
%% plotit = 1 would plot ft spectrum log scale (in dB), set = 0 for none
%% out:
%% filtSig2 = new filtered data (complex form)
%% datafft  = original fft
%% datafft2 = masked (convolved) fft

s=size(data);
% if exist('interpfac','var')~=1
%     interpfac = s(1);
% else
%     interpfac = s(1)*interpfac;
% end
if exist('plotit','var')~=1
    plotit=0;
end

% interpFac=max(s2,interpFac);
NFTL = max(2048,2^nextpow2(s(1)));
if filterCuts~=-1
    loW1=filterCuts(1);hiW1=filterCuts(2);
    loW2=filterCuts(3);hiW2=filterCuts(4);
    filtExp=1.0;  % 1 is pure hamming, <1 is broader, >1 is narrower;
    %     w=0.1; %weight for filter
    
    %%tmp=data(fRng,:)-repmat(mean(data),length(fRng),1); % Zero Mean
    %%s2=size(tmp);
    
    datafft = fft(data,NFTL,1);
    s2      = size(datafft);
    dataLength=s(1);
    %datafft=fft(data,[],1);
    
    nPtsf1=s2(1);
    if plotit==1
        figure(1111);
        tmp3=mean(abs(datafft),2);
        %         freqList = linspace(0,freq/2,floor(s2(1)/2)+1);
        %         plot(freqList,20*log10(tmp3((1:floor(s2(1)/2)+1))/max(tmp3(1:floor(s2(1)/2)+1)))); hold on;
        plot(linspace(-freq/2,freq/2,length(tmp3)),20*log10(fftshift(tmp3)+1e-9)); hold on;
        %         plot(linspace(-freq/2,freq/2,length(tmp3)),fftshift(20*log10(tmp3/max(tmp3(:))))); hold on;
    end
    
    loW1=max(round(nPtsf1*loW1/freq),1);
    hiW1=max(round(nPtsf1*hiW1/freq),1);
    loW2=max(round(nPtsf1*loW2/freq),2);
    hiW2=max(round(nPtsf1*hiW2/freq),2);
    mask2=zeros(nPtsf1,1);
    
    if loW1==hiW1
        mask2(loW1:loW2)=1;
    else
        %a=hamming(2*(hiW1-loW1)+1);
        a=hanning(2*(hiW1-loW1)+1);
        mask2(loW1:hiW1+1)=a(1:hiW1-loW1+2);
        mask2(hiW1+1:loW2)=ones(loW2-hiW1,1);
    end
    
    if loW2==hiW2
        mask2(loW2:hiW2)=1;
    else
        %a=hamming(2*(hiW2-loW2)+1).^filtExp;
        a=hanning(2*(hiW2-loW2)+1).^filtExp;
        mask2(loW2:hiW2)=a(hiW2-loW2+1:end);
    end
    if s2(2)>1
        mask2D=repmat(mask2,1,s2(2));
    else
        mask2D=mask2;
    end
    
    %
    %     mag         = masky3';
    %     phase       = angle(datafft);
    %     sig         = ifft(mag.*cos(phase)+ 1i*mag.*sin(phase));
    
    
    datafft2=datafft.*mask2D;
    %     datafft2=datafft + mask2D*2;
else
    if exist('fRng','var')~=1
        fRng=1:s(1);
    elseif length(fRng)<2 || fRng(2)<=fRng(1)
        fRng=1:s(1);
    end
    dataLength=length(fRng);
    NFTL = max(2048,2^nextpow2(dataLength));
    tmp=data(fRng,:)-repmat(mean(data),length(fRng),1); % Zero Mean
    s2=size(tmp);
    datafft=fft(tmp,NFTL,1);
    nPtsf1=s2(1);
    datafft(floor(nPtsf1/2):end,:)=0;
    datafft2 = datafft;
end
%filtSig2  = ifft(datafft2,[],1);
filtSig2  = ifft(datafft2,NFTL,1);
freqArray = linspace(0,freq/2,ceil(size(datafft2,1)/2));
filtSig2 = filtSig2(1:dataLength,:);

if plotit==1
    figure(1112);
    tmp3         = smooth(mean(abs(datafft2),2),9);
    freqList     = linspace(0,freq/2,floor(length(tmp3)/2)+1);
    plot(freqList,20*log10(tmp3(1:floor(length(tmp3)/2)+1)+1e-9)); hold on;
    
    [a,b]=max(tmp3);
    pkFq = freqList(b(1));
    %   plot(freqList,20*log10(tmp3((1:floor(s2(1)/2)+1))/max(tmp3(1:floor(s2(1)/2)+1)))); hold on;
end
end
