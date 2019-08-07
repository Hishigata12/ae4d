function HFData = TwoDNotch(HFData,Chan,Parms,Harms,FilterCutsSlow,FilterCutsFast)
%% TwoDNotch.m 
% Arguments-
%       HFData (MxN Array of Doubles): 1D-filtered HF Data
%       Chan (double): Channel Number to Analyze
%       Parms (struc): Structure of all parameters included in allparams/bScanParm
%       Harms (double): Number of Harmonics to filter (only takes out harmonics greater than what is set in FilterCuts
%       FilterCutsSlow (1x4 Vector of Doubles): Hanning Window Cutoffs in slow time (in Hz)
%       FilterCutsFast (1x4 Vector of Doubles): Hanning Window Cutoffs in fast time (in MHz)
% Output -
%       HFData (MxN Array of Doubles): 2D-filtered HF Data
%
% Written by Alex Alvarez 04/11/19; Updated 07/03/19
%% Plot pre-filtered 2D FFT
filtExp=1.0;  %  1 is pure hamming, <1 is broader, >1 is narrower;
s=size(HFData);
NFTL = max(2048,2^nextpow2(s(1)));
NFTL2 = max(2048,2^nextpow2(s(2))); 
datafft2 = fft2(real(HFData),NFTL,NFTL2); 
s2      = size(datafft2); 
figure(1111*Chan); 
tmp3=abs(datafft2); 
tmp4=20*log10(tmp3+1e-9);tmp4 = tmp4 - max(tmp4(:)); 
freq = Parms.daq.HFdaq.pulseRepRate_Hz; 
freq2 = Parms.daq.HFdaq.fs_MHz; 
x_ax = linspace(-freq/2,freq/2,s2(2));
y_ax = linspace(-freq2/2,freq2/2,s2(1)); 
subplot(2,1,1) 
imagesc(x_ax,y_ax,fftshift(tmp4),[-80 0]); colormap('gray'); xlabel('Slow Time Freq (Hz)'); ylabel('Fast Time Freq (MHz)'); colorbar; % AMA Addition 190410

SlowPts=s2(2); 
FastPts=s2(1); 
%% Create 2D Masks for each Harmonic given by FilterCutsSlow and FilterCutsFast
for p = 1:Harms
    lowslow1 = p*FilterCutsSlow(1); lowslow2 = p*FilterCutsSlow(2); highslow1 = p*FilterCutsSlow(3); highslow2 = p*FilterCutsSlow(4); % AMA Addition 190416 
    lowfast1 = p*FilterCutsFast(1); lowfast2 = p*FilterCutsFast(2); highfast1 = p*FilterCutsFast(3); highfast2 = p*FilterCutsFast(4); % AMA Addition 190416
    
    lowslow1=round(SlowPts*lowslow1/freq); 
    lowslow2=round(SlowPts*lowslow2/freq);
    highslow1=round(SlowPts*highslow1/freq); 
    highslow2=round(SlowPts*highslow2/freq);      
    
    lowfast1=round(FastPts*lowfast1/freq2); 
    lowfast2=round(FastPts*lowfast2/freq2); 
    highfast1=round(FastPts*highfast1/freq2); 
    highfast2=round(FastPts*highfast2/freq2);               
    
    mask2 = ones(FastPts,SlowPts);
    a=hanning(2*(lowslow2-lowslow1)+1); a =1-a;
    b=hanning(2*(highslow2-highslow1)+1).^filtExp; b = 1-b;
    c=hanning(2*(lowfast2-lowfast1)+1); c = 1-c; 
    d=hanning(2*(highfast2-highfast1)+1).^filtExp; d=1-d;
    SlowTmp = numel(a(1:lowslow2-lowslow1+1))+numel(zeros(highslow1-lowslow2-1,1))+numel(b(highslow2-highslow1+1:end));
    FastTmp = numel(c(1:lowfast2-lowfast1+1))+numel(zeros(highfast1-lowfast2-1,1))+numel(d(highfast2-highfast1+1:end));
    tmpa = repmat(a(1:lowslow2-lowslow1+1),1,FastTmp)'; 
    tmp1 = repmat(zeros(highslow1-lowslow2-1,1),1,FastTmp)'; 
    tmpb = repmat(b(highslow2-highslow1+1:end),1,FastTmp)'; 
    tmpc = repmat(c(1:lowfast2-lowfast1+1),1,SlowTmp); 
    tmp2 = repmat(zeros(highfast1-lowfast2-1,1),1,SlowTmp);          
    tmpd = repmat(d(highfast2-highfast1+1:end),1,SlowTmp); 

    tmpab = cat(2,tmpa,tmp1,tmpb);     
    tmpcd = cat(1,tmpc,tmp2,tmpd);  

    tmp = tmpab.*tmpcd;      

    mask2(lowfast1:highfast2,lowslow1:highslow2) = tmp; 
    mask2(end-highfast2:end-lowfast1,end-highslow2:end-lowslow1) = tmp;
    mask2D = mask2; 
    datafft3=datafft2.*mask2D;
    datafft2 = datafft3; 
    clear mask2; clear mask2D;
end 
%% Plot 2D-filtered FFT
tmp5=abs(datafft3); 
tmp6=20*log10(tmp5 + 1e-9);tmp6 = tmp6 - max(tmp6(:));
figure(1111); subplot(2,1,2)
imagesc(x_ax,y_ax,fftshift(tmp6),[-80 0]); colormap('gray'); xlabel('Slow Time Freq (Hz)'); ylabel('Fast Time Freq (MHz)'); colorbar; % AMA Addition 190410
%% IFFT for Filtered Data
TwoD_FiltData = ifft2(datafft3,NFTL,NFTL2); 
TwoD_FiltData = TwoD_FiltData(1:s(1),1:s(2));

HFData = TwoD_FiltData; 

end