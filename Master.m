[file, path] = uigetfile(fullfile(pwd,'*_info.dat'));
global floc
global nPoints
global m_burst_pts
floc = [path file];
[daq velmex hp] = call_file(1);
nPoints = velmex.XNStep*velmex.YNStep;
m_burst_pts = daq.duration_ms*daq.burstRepRate_Hz/1000;

%**************************************************************************
%Gets out HF and LF signals

n = ceil(nPoints/2); %Number of points in n depends on nPoints
[param HF LF] = read_ucsdi_data(floc,n);
if size(HF,3) == 2
    Sig.AE = HF(:,:,2); Sig.PE = HF(:,:,1);
else
    Sig.AE = HF;
end
Sig.AE_slow = Sig.AE';
Sig.LF = LF*50;
%**************************************************************************
%Creates original axes for current, x,y and z

Ax.m = linspace(0,daq.duration_ms,size(Sig.AE,2)); %Current axis
Ax.AE_dist_mm = size(Sig.AE,1)*1.48/daq.HFdaq.fs_MHz;
Ax.AE = linspace(0,Ax.AE_dist_mm,size(Sig.AE,1));
% Ax.PE_dist_mm = round(Ax.AE_dist_mm/2);
% Ax.PE = linspace(0,Ax.PE_dist_mm,size(Sig.PE,1));
Ax.x = linspace(-1*velmex.XDist,velmex.XDist,velmex.XNStep);
Ax.y = linspace(-1*velmex.YDist,velmex.YDist,velmex.YNStep);

%**************************************************************************
% Creates signals with different filter specifications
[Sig.filt.t.AE_none, Sig.filt.k.AE_none] = fft_and_filt(Sig.AE,daq.HFdaq.fs_MHz,'none');
[Sig.filt.t.AE_Hann, Sig.filt.k.AE_Hann, Ax.fft.F_axis, Sig.filt.k.filt_win] = fft_and_filt(Sig.AE,daq.HFdaq.fs_MHz,'Hann',[.5 .8 1.2 1.5]);
[Sig.filt.t.AE_Hann_low, Sig.filt.k.AE_Hann_low, Ax.fft.slow] = fft_and_filt(Sig.AE_slow,daq.HFdaq.pulseRepRate_Hz,'Hann',[50 75 850 950]);
[Sig.filt.t.AE_Hann_full, Sig.filt.k.AE_Hann_full] = fft_and_filt(Sig.filt.t.AE_Hann',daq.HFdaq.pulseRepRate_Hz,'Hann',[50 75 850 950]);
Ax.fft.flen = length(Ax.fft.F_axis); %flen2 = length(f2);
%AEf = AE_f+AE_f2;
[Sig.filt.t.AE_Wien, Sig.filt.k.AE_Wien] = fft_and_filt(Sig.AE,daq.HFdaq.fs_MHz,'Wiener',1,1,0.2);
[Sig.filt.t.AE_Wien_full, Sig.filt.k.AE_Wien_full] = fft_and_filt(flipud(Sig.filt.t.AE_Wien'),daq.HFdaq.pulseRepRate_Hz,'Wiener',1,600,300);
%[Sig.filt.t.AE_Wien2, Sig.filt.k.AE_Wien2] = fft_and_filt(Sig.AE,daq.HFdaq.fs_MHz,'Wiener',23,1,0.2);
[Sig.filt.t.AE_Hann, Sig.filt.k.AE_Hann, Ax.fft.F_axis, Sig.filt.k.filt_win] = fft_and_filt(Sig.filt.t.AE_Hann_low,daq.HFdaq.fs_MHz,'Hann',[.5 .8 1.2 1.5]);

%**************************************************************************
%Experimental basebanding
% [Sig.filt.bb.t.AE_Hann, Sig.filt.bb.k.AE_Hann1, Sig.filt.bb.k.AE_Hann2, Ax.bb.t_ax, Ax.bb.k_ax] = baseband(Sig.filt.t.AE_Hann,1,20);
%
% [Sig.filt.bb.t.AE_Hann_f, Sig.filt.bb.k.AE_Hann1_f, Sig.filt.bb.k.AE_Hann2_f, Ax.bb.t_ax, Ax.bb.k_ax] = baseband(Sig.filt.t.AE_Hann,1,20);
[Sig.filt.bb.t.AE, Sig.filt.bb.k.AE, Ax.bb.t_ax, Ax.bb.k_ax] = baseband2(Sig.filt.t.AE_Hann,1,daq.HFdaq.fs_MHz); % This basebanding works

%[Sig.filt.bb.t.AE_first, Sig.filt.bb.k.AE_first] = guassfilt(Sig.filt.bb.k.AE1,Ax.bb.k_ax,0.3);
% Filter Low after Baseband
[Sig.filt.t.AE_full_bb, Sig.filt.k.AE_full_bb] = fft_and_filt(Sig.filt.bb.t.AE',daq.HFdaq.pulseRepRate_Hz,'Hann',[100 150 850 950]);



Sig.forB = fliplr(Sig.filt.t.AE_Hann_full);
Ax.m = linspace(1,13,64);
X = abs(Sig.filt.t.AE_Hann_full(:,22));
X2 = abs(Sig.filt.t.AE_Hann_full(:,62));
X3 = abs(Sig.filt.t.AE_Hann_full(:,42));
plot(Ax.AE,X); ylim([0,.004])
a.snr.sig = X(460:520); a.snr.noise = X(300:360);
a.snr.sigrms = rms(a.snr.sig); a.snr.noiserms = rms(a.snr.noise);
a.snr.db = 20*log10(a.snr.sigrms/a.snr.noiserms);
a.snr.db2 = 20*log10(a.snr.sigrms/0.00011);
a.snr.I = max(Sig.LF)-min(Sig.LF);
a.ae_peak = max(a.snr.sig);
a.snr.sig2 = X2(460:520);
a.snr.sig3 = X3(460:520);
a.snr.db2a = 20*log10(rms(a.snr.sig2)/0.00011);
a.snr.db2b = 20*log10(rms(a.snr.sig3)/0.00011);
a.ae_peak2 = max(a.snr.sig2);
a.ae_peak3 = max(a.snr.sig3);

%%
%**************************************************************************
%B mode stuff

%z1 = max(max(real(Sig.filt.t.AE_Hann_full(300:end-300,:))));
%z2 = find(real(Sig.filt.t.AE_Hann_full(300:end-300))>z1*0.999,1);
%z = mod(z2,size(Sig.AE,1));% This needs to be fixed for accuracy
z = find(Ax.AE>65,1);

[Sig.b.Vxt, Sig.b.Vdx, Sig.b.Vdz, Sig.b.Vdt, Sig.b.Vdxdz, Sig.b.Vdxdt, Sig.b.Vdzdt, Sig.b.Vdxdzdt] = One_D_bmode(floc,z,'X',1);

%test = reshape(Sig.b.Vxt,[size(Sig.b.Vxt,3),size(Sig.b.Vxt,1),size(Sig.b.Vxt,2)]); %Verify if this is the same as the below reshaping
% for i = 1:size(Sig.b.Vxt,1)
%     for j = 1:size(Sig.b.Vxt,2)
%         for k = 1:size(Sig.b.Vxt,3)
%             test(k,i,j) = Sig.b.Vxt(i,j,k);
%         end
%     end
% end
% for i = 1:size(Sig.b.Vdt,1)
%     for j = 1:size(Sig.b.Vdt,2)
%         for k = 1:size(Sig.b.Vdt,3)
%             test2(k,i,j) = Sig.b.Vdxdzdt(i,j,k);
%             test3(k,i,j) = Sig.b.Vdt(i,j,k);
%         end
%     end
% end
% Sig.b.Vdt = test3;
% Sig.b.Vdxdzdt = test2;
% Sig.b.Vxt = test;
% clear test*;
%[Xq,Yq,Zq] = meshgrid(min(:(1/(size(Sig.b.Vdxdt,1)*3)):1,0:(1/(size(Sig.b.Vdxdt,2)*3)):1,0:(1/(size(Sig.b.Vdxdt,3)*3)):1);
for i = 1:size(Sig.b.Vxt,2)
    a1 = Sig.b.Vxt(:,i,:);
    a2 = squeeze(a1);
    imagesc(a2)
end

Smooth = interp3(Sig.b.Vdxdt,('spline'));
SVxt = interp3(Sig.b.Vxt,'spline');
s1 = SVxt(500,:,:);
s2 = squeeze(s1);
imagesc(s2)
colormap(jet)
%%
%imagesc(Sig.b.Vdxdt); %Need to fix this

% Ax.b.x = interp1(Ax.x, size(Smooth,1));
% Ax.b.t = interp1(Ax.m);
% Ax.b.z = [z-150 z+150]/daq.HFdaq.fs_MHz*1.48;


%Sig.b.Vxt2 = Sig.b.Vxt(50:300,1:16,:);%Use to zoom in
B_mode_disp(Smooth,1:100,10);

Ax.db = mean(mean(mean(SVxt)));
s1 = SVxt(320,:,:);
s2 = squeeze(s1);
for i = 1:size(s2,1)
    for j = 1:size(s2,2)
        
            smooth_db(i,j) = db(s2(i,j)/Ax.db);
        
    end
end
imagesc(s2)
colormap(jet)


% Ax.db = mean(mean(mean(Sig.b.Vdxdt)));
% for i = 1:size(Smooth,1)
%     for j = 1:size(Smooth,2)
%         for k = 1:size(Smooth,3)
%             smooth_db(i,j,k) = db(Smooth(i,j,k)/Ax.db);
%         end
%     end
% end
% smooth_db(smooth_db < -10) = -10;
% smooth_db(smooth_db > 10) = 10;
% s1 = smooth_db(300,:,:);
% s2 = squeeze(s1);
% imagesc(s2)
% colormap(jet)
%

% caxis([0 max(max(smooth_db))]);
%
% if min(min(smooth_db)) == -Inf || max(max(smooth_db)) == Inf
% smooth_db = medfilt2(smooth_db,[3 3]);
% end


%% Displays images of data
%**************************************************************************
test = Sig.filt.t.AE_full_bb;
test2 = abs(test);
test3 = Sig.filt.t.AE_Hann_low;
test5 = Sig.filt.bb.t.AE_first;
test4 = fliplr(Sig.filt.t.AE_Hann_full);
test4abs = abs(test4);

Ax.filt.AE_axis = linspace(0,Ax.AE_dist_mm,length(Sig.filt.t.AE_none));
%AE_bb = baseband_russ3(AEf,20,0.9);

%surf(m_axis,AE_axis,real(AE_bb)+1); hold;
figure;
imagesc(Ax.m,Ax.AE,real(Sig.filt.t.AE_Hann_full));
center_axis(real(Sig.filt.t.AE_Hann_full),Ax.AE,'HF',[30 50],'lin');
figure;
imagesc(Ax.m,Ax.AE,real(Sig.filt.bb.t.AE));
center_axis(real(Sig.filt.bb.t.AE),Ax.AE,'HF',[30 50],'lin');
figure;
imagesc(Ax.m,Ax.AE,abs(real(Sig.filt.t.AE_Wien_full)));
center_axis(real(Sig.filt.t.AE_Wien_full),Ax.AE,'HF',[50 100],'lin');
colormap(jet)
% figure;
% imagesc(Ax.m,Ax.AE,real(Sig.filt.t.AE_none));
% center_axis(real(Sig.filt.t.AE_none),Ax.AE,'HF',[50 80],'lin');
% figure;
% imagesc(Ax.m,Ax.AE,real(aelf3));
% center_axis(real(aelf3),Ax.AE,'HF',[50 80],'lin');
% figure;
% imagesc(Ax.m,Ax.bb.t_ax,real(Sig.filt.bb.t.AE_Hann));
% center_axis(real(Sig.filt.bb.t.AE_Hann),Ax.bb.t_ax,'HF',[50 80],'lin');
figure;
imagesc(Ax.m,Ax.AE,real(Sig.forB));
center_axis(real(Sig.forB),Ax.AE,'HF',[60 90],'lin');

figure;
imagesc(Ax.m,Ax.AE,real(test4abs));
center_axis(real(test4abs),Ax.AE,'HF',[50 90],'lin');
colormap(jet)

figure;
imagesc(Ax.m,Ax.AE,real(test2));
center_axis(real(test2),Ax.AE,'HF',[60 90],'lin');

test6 = Sig.filt.t.AE_none;
test9 = real(test) + imag(test);

figure;
imagesc(Ax.m,Ax.AE,real(test9));
center_axis(real(test9),Ax.AE,'HF',[60 90],'lin');

figure;
imagesc(Ax.m,Ax.AE,real(Sig.AE));
center_axis(real(Sig.AE),Ax.AE,'HF',[60 90],'lin');


% Mean Filter
for i = 5:128
    for j = 1:4096
        mt(j,i-4) = (1/5)*(test3(j,i-4)+test3(j,i-3)+test3(j,i-2)+test3(j,i-1)+test3(j,i));
    end
end




