clear;
[file, path] = uigetfile(fullfile(pwd,'*_info.dat')); %Gets file location
param = read_ucsdi_info([path file]); %Gets scan parameters
[~,~,LF] = read_ucsdi_data([path file],1); %Gets input current waveform
[file2, path2] = uigetfile(fullfile(pwd,'*PEParm.mat')); %gets US pulse waveform
PE = open([path2 file2]);
US = PE.TW.Wvfm1Wy;
ax.HFfreq = linspace(0,param.daq.HFdaq.fs_MHz,param.daq.HFdaq.pts); %Creates fast frequency axis
ax.LFfreq = linspace(0,param.daq.HFdaq.pulseRepRate_Hz,param.daq.HFdaq.NoBurstTriggers); %creates slow frequency axis


%N = param.velmex.XNStep;
%**************************************************************************
%Builds 4D Matrix***************~~~~~~~~~~~~~~~~**************
[~, HF1] = full_signal([path file],param,2); %Gets the raw data

X = w_slow_filt2(param,HF1,LF,1,[100 300]); %Filters in slow time

X = w_ae_filt2(param,X,US,1,[2 4]); %Filters in fast time

for i = 1:size(X,1)
    for j = 1:size(X,2)
        for k = 1:size(X{1,1},2)
             X2{i,j}(:,k) = envelope(real(X{i,j}(:,k))); %Converts to envelope signal 
        end
    end
end
M.Depth = 1.48*param.daq.HFdaq.pts/param.daq.HFdaq.fs_MHz;
M.x = linspace(0,M.Depth,param.daq.HFdaq.pts);
M.xL = find(M.x < 50);
M.xH = find(M.x > 30);
M.xT = intersect(M.xL,M.xH);
ax.depth = linspace(0,M.Depth,param.daq.HFdaq.pts);
ax.x = linspace(-param.velmex.XDist/2,param.velmex.XDist/2,param.velmex.XNStep);
%t = find(M.x > 45,1);

% O = real(X{13,1}(:,32));
% plot(ax.HFfreq,abs(fft(O)))



%plot_B_mode(real(O),[10 20],[600 1400],21:40)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = linspace(0,size(X{1},1)/param.daq.HFdaq.fs_MHz,size(X{1},1));
%%%%%%%%%%%%%%%%%%%%%%%%
b = waitbar(0,'Basebanding');
% for j = 1:size(X,1)
%     for k = 1:size(X,2)
%         for m = 1:size(X{1},2)
%             X_bb{j,k}(:,m) = X{j,k}(:,m).*(sqrt(2).*exp(1i*2*pi*2.7*t))';
%         end
%     end
%     waitbar(j/size(X,1)/2)
% end
waitbar(0.5,b,'Enveloping and converting to 4D matrix')
for i = 1:param.velmex.XNStep
    for j = 1: param.velmex.YNStep
       % HF((i-1)*sL+j,:,:) = HF1{i,j}; %converts cells to pseudo 4-D array
       HF(i,j,:,:) = envelope(real(X{i,j})); %Converts cell array to double
      % HF_bb(i,j,:,:) = X_bb{i,j};
    end
    waitbar(.5+i/param.velmex.XNStep/2,b,'Finalizing 4D array construction');
end

XdB = real(20*log10(real(HF)./max(max(max(real(HF(:,:,M.xT(1):M.xT(end),:)))))));
%XdB_bb = real(20*log10(real(HF_bb)./max(max(max(real(HF_bb))))));
delete(b)

%vidwrite(real(XdB),ax.depth,ax.x,M.xT(1),M.xT(end),1,size(X,1),25,25,[-8 0]); %Creates standard video before iradon
%%
s = 20;
b = waitbar(0,'Computing inverse radon transform');
for i = 1:s
    q = permute(squeeze(HF(:,1,:,i)),[2 1]);
 %   q2 = cRadFilt(q);
%R(:,:,i) = cBackProj(q2,pi/4,2,100); % Performs back projection operation to original image
R(:,:,i) = iradon(q,linspace(-30,30,321),'None');
waitbar(i/s,b,'Computing inverse radon transform');
end
delete(b)
ax.radonx = linspace(-25,25,size(R,1)); %creates x and z axes
%%
XdB = real(20*log10(R./max(max(max(R)))));
J = permute(R,[2 1 3]);
L = J - min(min(min(R)));
I = L./max(max(max(L)));
figure; 
for i = 1:s
    %I2 = histeq(I(:,:,i));
    imshow(I(:,:,i),[0 0.1])
    drawnow
end
    
%%
vidwrite(XdB,ax.depth,ax.radonx,[200 1200],[200 1200],[1 s],[-9 0]) %plots output signal
%imagesc(R)

