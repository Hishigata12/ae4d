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
%Gets index for number of points to shift for each average
switch Wave
    case 'Pulse'
        Avg = size(LF,2);
        xInd = size(HF,1);
        yInd = size(HF,2);
        for i = 1:Avg  
            test = corr(Temp,LF(:,i));
            if test > 0.1
                cor = xcorr(Temp,LF(:,i));
                s = round(length(cor)/2)+find(cor == max(cor));
                if s~=0
                    LF_Shift(:,i) = circshift(LF(:,i),s);
                    HF_s = round(s*HFTimeSamp/LFSamp);
                    for j = 1:xInd
                        for k = 1:yInd
                            HF_Shift(j,k,:,:,i) = circshift(HF(j,k,:,:,i),HF_s,4);
                        end
                    end
                end
            end
            multiWaitbar('Shifting LF and HF in time',i/Avg);
        end
end

LF2 = LF_Shift;
HF2 = HF_Shift;
        
        




