function [LF2, HF2] = IndividualAvg(param, LF, HF, Temp, Wave)

%Description
%Aligns signal between multiple averages

%Inputs
%param is the paramater file associated with the AE scan
%LF is the low frequency matrix
%LF will be used as the main template for alignment
%HF is the full (X x Y x Depth x Slow Time x Avg) AE matrix given a
%particular channel
%Temp is the channel number for LF used as template
%Wave is the type of waveform for specific processing

%Gets universal parameters for HF and LF conversions
LFSamp = param.daq.LFdaq.fs_Hz/1e3;
HFDepthSamp = param.daq.HFdaq.fs_MHz*1e3;
HFTimeSamp = param.daq.HFdaq.pulseRepRate_Hz/1e3;
HF = permute(HF,[1 2 4 3 5]);
%Gets index for number of points to shift for each average
switch Wave
    case 'Pulse'
        Avg = size(LF,2);
        xInd = size(HF,1);
        yInd = size(HF,2);
        for j = 1:xInd
            for k = 1:yInd
                for i = 1:Avg
                    test = corr(Temp,LF(:,i));
                    if abs(test) < 0.8
                        Temp = (Temp - mean(Temp))/std(Temp);
                        Temp2 = LF(:,i);
                        Temp2 = (Temp2 - mean(Temp2))/std(Temp2);
                        cor = xcorr(Temp,Temp2);
                        s = round(length(cor)/2)+find(abs(cor) == max(abs(cor)));
                        if s~=0
                            LF_Shift(:,i,j,k) = circshift(LF(:,i),s);
                            HF_s = round(s*HFTimeSamp/LFSamp);
                            %                     figure(i);
                            %                     plot(LF_Shift(:,i));
                            
                            HF_Shift(j,k,:,:,i) = circshift(HF(j,k,:,:,i),HF_s,4);
                        end
                    else
                        LF_Shift(:,i) = LF(:,i);
                        HF_Shift(j,k,:,:,i) = HF(j,k,:,:,i);
                    end
                end
            end
            multiWaitbar('Shifting LF and HF in time',j/xInd);
        end
    case 'ECG'
        Avg = size(LF,2);
        xInd = size(HF,1);
        yInd = size(HF,2);
        b = (1/51)*ones(51,1);
        a = 1;
        for j = 1:xInd
            for k = 1:yInd
                c = 1;
                for i = 1:Avg
                    %                     LF_filt = filtfilt(b,a,LF(:,i)); %For smoothing the Sig
                    LF_filt = LF(:,i);
                    peaks = findpeaks(LF_filt);
                    s = find(LF_filt == max(peaks),1);
                    temp_max = abs(max(Temp));
                    
                    s_temp = find(abs(Temp) == temp_max);
                    shift = s_temp - s;
                        T2 = circshift(LF(:,i),shift);
                        S2 = find(abs(T2) == max(abs(T2)));
                        if abs(S2-s_temp) < 3
                            LF_Shift(:,c) = T2;
                            
                            HF_s = round(shift*HFTimeSamp/LFSamp);
                            %                     figure(i);
                            %                     plot(LF_Shift(:,i));
                            
                            HF_Shift(j,k,:,:,c) = circshift(HF(j,k,:,:,i),HF_s,4);
                        end
                        c = c+1;
                end
            end
            %             end
            multiWaitbar('Shifting LF and HF in time',i/Avg);
        end
end
switch ndims(LF_Shift)
    case 4
        LF_Shift = mean(LF_Shift,4);
        LF_Shift = mean(LF_Shift,3);
    case 3
        LF_Shift = mean(LF_Shift,3);
end

figure(7); plot(LF_Shift);

HFT = squeeze(HF_Shift(1,1,floor(size(HF,3)/2),:,1));
HFT = (HFT - mean(HFT))/std(HFT);
for i = 1:size(HF,1)
    HFC = squeeze(HF_Shift(i,1,floor(size(HF,3)/2),:,1));
    HFC = (HFC - mean(HFC))/std(HFC);
    HFC2 = squeeze(HF_Shift(i,1,floor(size(HF,3)/1.5),:,1));
    HFC2 = (HFC2 - mean(HFC2))/std(HFC2);
    HFC3 = squeeze(HF_Shift(i,1,floor(size(HF,3)/3),:,1));
    HFC3 = (HFC3 - mean(HFC3))/std(HFC3);
    HFC4 = squeeze(HF_Shift(i,1,floor(size(HF,3)/4),:,1));
    HFC4 = (HFC4 - mean(HFC4))/std(HFC4);
    cor = xcorr(HFT,HFC);
    cor2 = xcorr(HFT,HFC2);
    cor3 = xcorr(HFT,HFC3);
    cor4 = xcorr(HFT,HFC4);
    s = round(length(cor)/2)+find(abs(cor) == max(abs(cor)));
     s2 = round(length(cor2)/2)+find(abs(cor2) == max(abs(cor2)));
      s3 = round(length(cor3)/2)+find(abs(cor3) == max(abs(cor3)));
       s4 = round(length(cor4)/2)+find(abs(cor4) == max(abs(cor4)));
       s5 = round(mean([s s2 s3 s4]));
    HF_Shift(i,:,:,:,:) = circshift(HF_Shift(i,:,:,:,:),s5,4);
end


LF2 = mean(LF_Shift,2);
HF2 = mean(HF_Shift,5);
        
        
% figure;
% for i = 1:size(HF,1)
%      imagesc(squeeze(HF_Shift(i,1,:,:,1)));
% figure;
% plot(squeeze(LF_Shift(:,i,:)));
%      drawnow;
%      pause(0.2)
% end


