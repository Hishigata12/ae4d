function [HF,LF,PE,LF1,Location,NumThrown,NumSplit] = PostAverageZeta_v2(Parms,LF1,Location,SingleAvg,HFSingle,filterCuts,ScanPt,XPt,YPt,Chan,chIdx,Manually,ShowPlots,DataToLoad,LoadPE,PEPRF,PESingle)
%% Give correlation coefficient after each average
%% Arguments
    % Parms (struc)- Parameter file (typically called from bScanParm or equivalent)
    % LF1 (Nx1 double or []) - Low frequency waveform to use for shifting
    % Location (double or []) - Sample that indicates center of waveform
    % SingleAvg (2D Array of Double) - Low frequency data set (SlowTime x Average)
    % HFSingle (3D Array of Doubles) - High frequency data set new (FastTime x SlowTime x Average)   
    % filtercuts (4x1 double) - [Low1 Low2 High1 High2] - Frequency cutoffs for Hanning window
    % ScanPt (double) - Scan Point
    % XPt (double) - X Scan Point
    % YPt (double) - Y Scan Point
    % Chan (double) - Channel Number
    % chIdx (double) - Channel Index
    % Manually (Boolean) - Flag for manually removing and splitting averages
    % ShowPlots (Boolean) - Flag for showing plots during cross-correlation process
    % DataToLoad (String) - Will signal to load 'UnShift,' 'UnAvg,' 'Avg,' or 'Filt' data
    % LoadPE (Boolean) - Flag for Loading PE
    % PEPRF (double) - the frame rate of PE
    % PESingle (4D Array of Doubles) -  PE data set (Fast Time x Slow Time(NumPulses) x Elements x Averages)
%% Outputs
    % HF (2-3D array of double) - HF Data Set new - Fast Time x Slow Time x Num of Avg (if DataToLoad is 'UnAvg')
    % LF (1-2D array of double) - LF Data Set - Amplitude x Num of Avg (if DataToLoad is 'UnAvg')
    % PE (3-4D array of double) - PE RF Data Set - Fast Time x Slow Time x Elements x Num of Avg (if DataToLoad is 'UnAvg)
    % LF1 (Nx1 double) - Low frequency waveform to use for shifting in remainder of scan points
    % Location (double) - Sample that indicates center of waveform
    % NumThrown (double) - Number of Averages thrown during the averaging
    % NumSplit (double) - Number of Averages split during the averaging
%% Initial definitions
AllAvgTime = tic; % Start a clock to see how long it takes
ConversionFactor = (Parms.Daq.HF.PulseRate*Parms.Scan.Duration_ms/1000)/Parms.Daq.LF.Samples; % Determine how many HF samples can be converted to LF Samples
if ~isempty(PEPRF)
    PEConversionFactor = (PEPRF*Parms.Scan.Duration_ms/1000)/Parms.Daq.LF.Samples;
end
if LoadPE && ~isempty(PESingle)
    PESingle = permute(PESingle,[2 1 3 4]);
end
HFSingle = permute(HFSingle,[2 1 3 4]);
LFChanns = Parms.Chan.LFChanns;
% Someday include something if HFChans ~= LFChans   
subplotrows = ceil(sqrt(Parms.Scan.Avg)); % Determine the number of rows needed for an appropriate subplot array for the number of averages
subplotcolumns = ceil(Parms.Scan.Avg/subplotrows); % Do the same with the columns
PeakWidth = 400;
NoiseCorrection = 1; % Correction for signals with large T/P waves
if ShowPlots % Create a figure if the flag for ShowPlots is set to 1
    OriginalPlots = figure(600); OriginalPlots.Units = 'normalized'; OriginalPlots.OuterPosition = [0 0 1 1]; clf; % Set the figure to take the entirety of the screen
end
for AvgNum = 1:Parms.Scan.Avg % Define the 
    LFSingle(:,AvgNum) = SingleAvg(:,AvgNum)-mean(SingleAvg(:,AvgNum)); % Remove any DC offset that might be present
    [pks,locs] = findpeaks(LFSingle(:,AvgNum),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(LFSingle(:,AvgNum)),'Annotate','extents');
    [negpks,neglocs] = findpeaks(-LFSingle(:,AvgNum),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(-LFSingle(:,AvgNum)),'Annotate','extents');
    if max(abs(pks)) >= max(abs(negpks))
        peaks = pks; locas = locs;
    else
        peaks = negpks; locas = neglocs;
    end
    if ShowPlots  % && numel(peaks) == 1
        subplot(subplotrows,subplotcolumns,AvgNum)
        plot(LFSingle(:,AvgNum)); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['Pt' num2str(ScanPt) ' Chan' num2str(LFChanns(Chan)) ' Avg' num2str(AvgNum)]); axis tight;  hold on
        plot(locs,pks, ' ro ');  
        plot(neglocs,-negpks,' bo '); hold off
    end
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
    CenterLag = location-round(Parms.Daq.LF.Samples/4); % CenterLag is now QuarterLag to include as many T's as possible
    SignalROI = location-PeakWidth:location+PeakWidth;
    if SignalROI(end) > size(LF1,1)
        SignalROI(end) = size(LF1,1)-1;
    end
    if SignalROI(1) < 1
        SignalROI(1) = 1;
    end
    NoiseLF = cat(1,LF1(1:SignalROI(1)),LF1(SignalROI(end)+1:end))./NoiseCorrection;
    MeanNoiseLF = mean(NoiseLF(:));
    StdNoiseLF = std(NoiseLF(:)); 
    if CenterLag < 0
        CenterLag = -CenterLag;
        PadLF1 = size(LF1,1)-size(LF1(1:end-CenterLag+1),1);        
        LFPad = MeanNoiseLF + StdNoiseLF.*randn(PadLF1,1);   
        LF1Shifted = LFPad;
        LF1Shifted(CenterLag:size(LF1,1)) = LF1(1:end-CenterLag+1);    
    elseif CenterLag > 0
        CenterLag = CenterLag;
        LF1Shifted = LF1(CenterLag+1:end,:);
        PadLF1 = size(LF1,1)-size(LF1Shifted,1);
        LFPad = MeanNoiseLF + StdNoiseLF.*randn(PadLF1,1);  
        LF1Shifted(size(LF1,1)-CenterLag+1:size(LF1,1),:) = LFPad; 
    else 
        LF1Shifted = LF1;
    end
    LF1 = LF1Shifted;
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
    Location = location;
else
    location = Location;
end
if ShowPlots
    figure(700); clf; figure(800); clf; % Open two new figures for plotting cross-correlations and shifted signals
end
LF = zeros(size(LF1)); % Sum the LF waveform
if size(HFSingle,3) >= 1
    HF = zeros(size(HFSingle(:,:,1))); % Sum the HF waveform 
else
    HF = zeros(size(HFSingle));
end
if size(PESingle,4) >= 1
    PE = zeros(size(PESingle(:,:,:,1)));
else
    PE = zeros(size(PESingle));
end
padLF = 40;
LF1 = LF1';
MeanLF1 = 'replicate';
LF1 = padarray(LF1',[padLF,0],MeanLF1,'both')';
LF1 = 2*real(fft_filt2D_new(LF1',[],Parms.Daq.LF.Rate_hz*1000,filterCuts,0))';
LF1 = LF1(:,1+padLF:end-padLF)';
for AvgNum = 1:NumAvgs
    LF2 = LFSingle(:,AvgNum); % Set the LF2 signal equal to each average (including the first)
    HF2 = HFSingle(:,:,AvgNum); % Create a dummy variable HF2 for shifting each average of data
    if LoadPE && ~isempty(PESingle)
        PE2 = PESingle(:,:,:,AvgNum);
    end
    LFFilt = LF2';
    MeanLFFilt = 'replicate';
    LFFilt = padarray(LFFilt',[padLF,0],MeanLFFilt,'both')';
    LFFilt = 2*real(fft_filt2D_new(LFFilt',[],Parms.Daq.LF.Rate_hz*1000,filterCuts,0))';
    LFFilt = LFFilt(:,1+padLF:end-padLF)';
    [cors, lags] = xcorr(abs(LF1),abs(LFFilt)); % Calculate the lags and correlation values between the LF1 and LF2 signal
    [~,I] = max(abs(cors)); % Find the index of the maximum absolute correlation
    lagDiff = lags(I);
    SignalROI = (location-lagDiff)-PeakWidth:(location-lagDiff)+PeakWidth;
    if SignalROI(end) > size(LF2,1)
        SignalROI(end) = size(LF2,1)-1;
    end
    if SignalROI(1) < 1
        SignalROI(1) = 1;
    end    
    NoiseLF = cat(1,LF2(1:SignalROI(1)),LF2(SignalROI(end)+1:end))./NoiseCorrection;
    MeanNoiseLF = mean(NoiseLF(:));
    StdNoiseLF = std(NoiseLF(:)); 
    
    ConvFactor = round(lagDiff*ConversionFactor); % Calculate how many HF samples this lag corresponds to    
    SignalROIHF = round(SignalROI.*ConversionFactor);
    NoiseHF = cat(1,HF2(1:SignalROIHF(1),:),HF2(SignalROIHF(end)+1:end,:))./NoiseCorrection;    
    MeanNoiseHF = mean(NoiseHF,1);
    StdNoiseHF = std(NoiseHF,0,1);

    if LoadPE && ~isempty(PESingle)
        PEConvFactor = round(lagDiff*PEConversionFactor);
        SignalROIPE = round(SignalROI.*PEConversionFactor);
        NoisePE = cat(1,PE2(1:SignalROIPE(1),:,:),PE2(SignalROIPE(end)+1:end,:,:))./NoiseCorrection; %% Might need to change what section the noise is taken from
        MeanNoisePE = mean(NoisePE,1);
        StdNoisePE = std(NoisePE,0,1);
    end
    %% Plot the correlations
    if ShowPlots
        figure(700);
        subplot(subplotrows2,subplotcolumns2,AvgNum)
        plot(lags,cors); xlabel('Lags'); ylabel('Correlation Value'); title(['Correlation between 1&' num2str(AvgNum)]); axis tight;      
        hold on
    end
    %% Shift Low and High Frequency Peaks   
    if lagDiff > 0
        PadLF2 = size(LF2,1)-size(LF2(1:end-lagDiff+1),1);
        LFPad = MeanNoiseLF + StdNoiseLF.*randn(PadLF2,1);  
        LFShifted = LFPad;
        LFShifted(lagDiff:size(LF2,1)) = LF2(1:end-lagDiff+1);

        PadHF2(1) = size(HF2,1)-size(HF2(1:end-ConvFactor+1,:),1);
        PadHF2(2) = size(HF2,2);
        HFPad = MeanNoiseHF + StdNoiseHF.*randn(PadHF2(1),PadHF2(2));
        HFShifted = HFPad;
        HFShifted(ConvFactor:size(HF2,1),1:size(HF2,2)) = HF2(1:end-ConvFactor+1,:);
        if LoadPE && ~isempty(PESingle)
            PadPE2(1) = size(PE2,1)-size(PE2(1:end-PEConvFactor+1,:,:),1);
            PadPE2(2) = size(PE2,2);
            PadPE2(3) = size(PE2,3);
            PEPad = MeanNoisePE + StdNoisePE.*randn(PadPE2(1),PadPE2(2),PadPE2(3));
            PEShifted = PEPad;
            PEShifted(PEConvFactor:size(PE2,1),1:size(PE2,2),1:size(PE2,3)) = PE2(1:end-PEConvFactor+1,:,:);
        end
    elseif lagDiff < 0
        lagDiff = -lagDiff;
        LFShifted = LF2(lagDiff+1:end,:);
        PadLF2 = size(LF2,1)-size(LFShifted,1);
        LFPad = MeanNoiseLF + StdNoiseLF.*randn(PadLF2,1); 
        LFShifted(size(LF2,1)-lagDiff+1:size(LF2,1),:) = LFPad; 

        ConvFactor = -ConvFactor;
        HFShifted = HF2(ConvFactor+1:end,:);
        PadHF2(1) = size(HF2,1)-size(HFShifted,1);
        PadHF2(2) = size(HF2,2);
        HFPad = MeanNoiseHF + StdNoiseHF.*randn(PadHF2(1),PadHF2(2));
        HFShifted(size(HF2,1)-ConvFactor+1:size(HF2,1),:) = HFPad;

        if LoadPE && ~isempty(PESingle)
            PEConvFactor = -PEConvFactor;
            PEShifted = PE2(PEConvFactor+1:end,:,:);
            PadPE2(1) = size(PE2,1)-size(PEShifted,1);
            PadPE2(2) = size(PE2,2);
            PadPE2(3) = size(PE2,3);
            PEPad = MeanNoisePE + StdNoisePE.*randn(PadPE2(1),PadPE2(2),PadPE2(3));
            PEShifted(size(PE2,1)-PEConvFactor+1:size(PE2,1),:,:) = PEPad;
        end
    else 
        LFShifted = LF2;
        HFShifted = HF2;
        if LoadPE && ~isempty(PESingle)
            PEShifted = PE2;
        end
    end
                           
    % Plot newly shifted LF waveforms
    if ShowPlots
        figure(800);
        subplot(subplotrows2,subplotcolumns2,AvgNum)
        plot(LFShifted); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['Pt' num2str(ScanPt) ' Chan' num2str(LFChanns(Chan)) ' Avg' num2str(AvgNum)]); axis tight;   
    end
    if strcmpi(DataToLoad,'Avg')
        tmp2 = '/LF'; % Create a temporary variable for an h5 saving variable
        tmp3 = '/HF'; % Repeat with HF
        tmp4 = [pwd, '\', Parms.Format.filename2, '_Raw4dUnAveraged.h5'];
        tmp6 = '/PE'; % Repeat with PE
        startLFUnAvg = [1 AvgNum chIdx XPt YPt];
        startHFUnAvg = [1 1 AvgNum chIdx XPt YPt];
        startPEUnAvg  = [1 1 1 AvgNum XPt YPt];
        countLFUnAvg = [Parms.Daq.LF.Samples,1,1,1,1];
        countHFUnAvg = [round(Parms.daq.HFdaq.pts),Parms.daq.HFdaq.NoBurstTriggers,1,1,1,1];
        countPEUnAvg = [size(PEShifted,2), size(PEShifted,1),size(PEShifted,3),1,1,1];
        HFShifted = permute(HFShifted, [2 1 3]);
        SaveUnAvgTime = tic;
        if LoadPE
            PEShifted = permute(PEShifted, [2 1 3 4]);
            h5write(tmp4,tmp6,PEShifted,startPEUnAvg,countPEUnAvg);
        end        
        h5write(tmp4,tmp2,LFShifted,startLFUnAvg,countLFUnAvg)
        h5write(tmp4,tmp3,HFShifted,startHFUnAvg,countHFUnAvg)
        TotalSaveUnAvgTime = toc(SaveUnAvgTime);
        disp([num2str(TotalSaveUnAvgTime),' Seconds to Save UnAveraged Data for chIdx=' num2str(chIdx) ' @ X =' num2str(XPt) ' Y=' num2str(YPt)])
        
        
        HFShifted = permute(HFShifted, [2 1 3]);
        if LoadPE
            PEShifted = permute(PEShifted, [2 1 3 4]);
        end  
        LF(:,AvgNum) = LFShifted; % Sum the LF waveform
        HF(:,:,AvgNum) = HFShifted; % Sum the HF waveform
        if LoadPE && ~isempty(PE)
            PE(:,:,:,AvgNum) = PEShifted; % Sum the PE image
        end
        [pks,locs] = findpeaks(LF(:,AvgNum),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(LF(:,AvgNum)),'Annotate','extents');
        [negpks,neglocs] = findpeaks(-LF(:,AvgNum),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(-LF(:,AvgNum)),'Annotate','extents');
        if max(abs(pks)) >= max(abs(negpks))
            peaks = pks; locas = locs;
        else
            peaks = negpks; locas = neglocs;
        end
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
            plot(LF(:,AvgNum)); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['Pt' num2str(ScanPt) ' Chan' num2str(LFChanns(Chan)) ' Avg' num2str(AvgNum)]); axis tight;  hold on
            plot(locs,pks, ' ro ');  
            plot(neglocs,-negpks,' bo '); hold off
        end                
    else
        LF(:,AvgNum) = LFShifted;
        HF(:,:,AvgNum) = HFShifted;
        if LoadPE && ~isempty(PE)
            PE(:,:,:,AvgNum) = PEShifted;
        end
    end
end

if strcmpi(DataToLoad,'Avg')
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
            [pks,locs] = findpeaks(LF(:,SplitAvg),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(LF(:,SplitAvg)),'Annotate','extents');
            [negpks,neglocs] = findpeaks(-LF(:,SplitAvg),'MinPeakDistance',0.1*Parms.Daq.LF.Samples,'MinPeakHeight',0.6*max(-LF(:,SplitAvg)),'Annotate','extents');        
            if (max(abs(pks)) >= max(abs(negpks)))
                peaks = pks; locas = locs;
            elseif isempty(negpks)
                peaks = pks; locas = locs;
            else
                peaks = negpks; locas = neglocs;
            end
            
            LF(:,end+1) = LF(:,SplitAvg);
            SignalROI(1,:) = locas(1)-PeakWidth:locas(1)+PeakWidth;
            SignalROI(2,:) = locas(2)-PeakWidth:locas(2)+PeakWidth;            
            SplitPoint = round((locas(2)-locas(1))/2); 
            FirstHalf = LF(1:SplitPoint,SplitAvg); 
            SecondHalf = LF(SplitPoint+1:end,SplitAvg); 
            NoiseLF = cat(1,LF(1:SignalROI(1,1),SplitAvg),LF(SignalROI(1,end)+1:SignalROI(2,1),SplitAvg),LF(SignalROI(2,end)+1:end,SplitAvg))./NoiseCorrection;
            MeanNoiseLF = mean(NoiseLF(:)); StdNoiseLF = std(NoiseLF(:)); 
            PadFirstHalf = size(LF(:,SplitAvg),1)-size(FirstHalf,1);
            FirstHalfPad = MeanNoiseLF + StdNoiseLF.*randn(PadFirstHalf,1);            
            LF(SplitPoint+1:end,SplitAvg)= FirstHalfPad;
            PadSecondHalf = size(LF(:,SplitAvg),1)-size(SecondHalf,1);
            SecondHalfPad = MeanNoiseLF + StdNoiseLF.*randn(PadSecondHalf,1);            
            LF(1:SplitPoint,end) = SecondHalfPad;            

            HF(:,:,end+1) = HF(:,:,SplitAvg);
            SignalROIHF = round(SignalROI.*ConversionFactor);
            SplitPointHF = round(SplitPoint.*ConversionFactor);
            FirstHalfHF = HF(1:SplitPointHF,:,SplitAvg);
            SecondHalfHF = HF(SplitPointHF+1:end,:,SplitAvg);            
            NoiseHF = cat(1,HF(1:SignalROIHF(1,1),:,SplitAvg),HF(SignalROIHF(1,end):SignalROIHF(2,1),:,SplitAvg),HF(SignalROIHF(2,end):end,:,SplitAvg))./NoiseCorrection;                        
            MeanNoiseHF = mean(NoiseHF,1); StdNoiseHF = std(NoiseHF,0,1);
            PadFirstHalfHF = size(HF(:,:,SplitAvg),1)-size(FirstHalfHF,1);
            FirstHalfHFPad = MeanNoiseHF + StdNoiseHF.*randn(PadFirstHalfHF,size(HF(:,:,SplitAvg),2),1);
            HF(SplitPointHF+1:end,:,SplitAvg) = FirstHalfHFPad;
            PadSecondHalfHF = size(HF(:,:,SplitAvg),1)-size(SecondHalfHF,1);
            SecondHalfHFPad = MeanNoiseHF + StdNoiseHF.*randn(PadSecondHalfHF,size(HF(:,:,SplitAvg),2),1);            
            HF(1:SplitPointHF,:,end) = SecondHalfHFPad;
            if LoadPE && ~isempty(PE)
                PE(:,:,:,end+1) = PE(:,:,:,SplitAvg);
                SignalROIPE = round(SignalROI.*PEConversionFactor);
                SplitPointPE = round(SplitPoint.*PEConversionFactor);
                FirstHalfPE = PE(1:SplitPointPE,:,:,SplitAvg);
                SecondHalfPE = PE(SplitPointPE+1:end,:,:,SplitAvg);
                PadFirstHalfPE = size(PE(:,:,:,SplitAvg),1)-size(FirstHalfPE,1);
                NoisePE = cat(1,PE(1:SignalROIPE(1,1),:,:,SplitAvg),PE(SignalROIPE(1,end):SignalROIPE(2,1),:,:,SplitAvg),PE(SignalROIPE(2,end):end,:,:,SplitAvg))./NoiseCorrection;
                MeanNoisePE = mean(NoisePE,1); StdNoisePE = std(NoisePE,0,1); 
                FirstHalfPEPad = MeanNoisePE + StdNoisePE.*randn(PadFirstHalfPE,size(PE(:,:,:,SplitAvg),2),1);
                PE(SplitPointPE+1:end,:,:,SplitAvg) = FirstHalfPEPad;
                PadSecondHalfPE = size(PE(:,:,:,SplitAvg),1)-size(SecondHalfPE,1);
                SecondHalfPEPad = MeanNoisePE + StdNoisePE.*randn(PadSecondHalfPE,size(PE(:,:,SplitAvg),2),1);            
                PE(1:SplitPointPE,:,:,end) = SecondHalfPEPad;     
            end
        end
    else
        NumSplit = 0;
    end
    if AvgExist ~= 0
        NumThrown = numel(AvgsToThrowAway);
        LF(:,AvgsToThrowAway) = []; LF = squeeze(LF);  % Set Averages that the user does not want to include in averaging to null and recreate size of arrzy
        HF(:,:,AvgsToThrowAway) = []; HF = squeeze(HF); % Repeat for the HF
        if LoadPE && ~isempty(PE)
            PE(:,:,:,AvgsToThrowAway) = []; PE = squeeze(PE);
        end
    else
        NumThrown = 0;
    end

    LF = sum(LF,2);
    HF = sum(HF,3);
    PE = sum(PE,4);
    LF = LF./NumAvgs; % Average the LF waveform
    HF = HF./NumAvgs; % Average the HF waveform
    HF = permute(HF, [2 1]);
    if LoadPE && ~isempty(PE)
        PE = PE./NumAvgs; % Average the PE image
        PE = permute(PE, [2 1 3]);
    else
        PE = [];
    end
    figure(900+Chan);
    plot(LF); xlabel('LF Samples'); ylabel('Amplitude (V)'); title(['Averaged  ' ' LF Chann ' num2str(LFChanns(Chan))]); axis tight; % Eventually will be LFAveraged(:,XPt,YPt,Chan)     
    hold on
    disp([num2str(toc(AllAvgTime)-TotalSaveUnAvgTime),' Seconds to Average Scan for chIdx=' num2str(chIdx) ' @ X=' num2str(XPt) ' Y=' num2str(YPt)])
else
    HF = permute(HF, [2 1 3]);
    if LoadPE && ~isempty(PESingle)
        PE = permute(PE, [2 1 3 4]);
    else
        PE = [];
    end
end
end