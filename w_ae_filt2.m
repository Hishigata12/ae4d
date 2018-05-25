%x contains output filtered signal
%HF contains 4D AE data
%US contains the US waveform
%n is the point along slow time that contains the noise spectrum in the
%depth direction

function y = w_ae_filt2(param,HF,us,mode,c)

%mode = 0;
%~~~~Create Filter~~~~%
FsUS = 250;
FsAE = param.daq.HFdaq.fs_MHz;
Lus = length(us);
Lae = param.daq.HFdaq.pts;
if Lus > Lae
    us = us(1:Lae);
    Lus = Lae;
end




frequs = linspace(0,FsUS,Lus);
freqae = linspace(0,FsAE,Lae);
max_fUS = length(frequs)/2;
max_fAE = length(freqae)/2;
fus = frequs(1:end/2);
fae = freqae(1:end/2);
% q = find(fus>15,1); % Fix US samples to AE
% US = fft(us);
% USs = US(1:end/2);
% US4 = US(1:q);
% US5 = interp1(1:q,US4,linspace(0,q,max_fUS));


% Hs = H(1:end/2);
% d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1.5/FsAE*2,2/FsAE*2,3.5/FsAE*2,4/FsAE*2,60,1,60);
% Hd = design(d,'equiripple');


%~~~~~Filter in fast time~~~~%
fspec = '* %d percent complete \n*';
HF_xy = size(HF);
HF_zt = size(HF{1});

%%%%% FREQUENCY FILTERING %%%%%%%%%

if mode == 0
    
    Fc1 = find(fae > c(1),1);
    Fc2 = find(fae > c(2),1);
    for i = 1:Lae
        if (i >= Fc1 && i <= Fc2) || (i >= Lae-Fc2 && i <= Lae-Fc1)
            Hd(i) = 1;
        else
            Hd(i) = 0;
        end
    end
    
    hd =ifft(Hd);
    w1q = hamming(100);
    wq = circshift(w1q,100/2);
    h = padarray(wq,size(HF{1},1)-100,'post').*hd';
    H = fft(h,size(HF{1},1));
    
        
elseif mode == 1
    if FsUS > FsAE
        RefPulse       = resample(us,FsAE,FsUS);
   
    else
        RefPulse = resample(us,FsUS,FsAE);
    end
    
    RefPulse = RefPulse/(sum(abs(RefPulse)));
    H = fft(RefPulse);
    H = interp1(linspace(0,FsAE,length(RefPulse)),H',linspace(0,FsAE,Lae));
    a = find(freqae > 0.4,1);
    s = length(H);
    H2(1:(a-1)) = 0;
    H2(a:s-(a-1)) = 1;
    H2(s-(a-1):s) = 0;
    H = H.*H2;
end 
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            % X{i,j} = HF{i,j}(1:end/2,:);
            X2{i,j} = zeros(HF_zt(1),HF_zt(2));
            y{i,j} = zeros(HF_zt);
        end
    end
    
   % fprintf('Filtering 4D data\n')
    b = waitbar(0,'Filtering 4D data');
    if size(H,2) > size(H,1)
        H = H';
    end
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            for k = 1:HF_zt(2)
                X2{i,j}(:,k) = fft(HF{i,j}(:,k)).*H;
                %X2(i,j,k,:) = filter(Hd,squeeze(X(i,j,k,:)));
            end
        end
        waitbar(i/HF_xy(1)/2,b,'Fast Time Filtering')
    end
    
    
    %fprintf('Filtering 4D data\n')
    
    %delete(b)
    % Xf = flip(X2,3);
    % Xs = cat(3,X2,Xf);
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            for k = 1:HF_zt(2)
                y{i,j}(:,k) = ifft(X2{i,j}(:,k),HF_zt(1));
            end
        end
        waitbar(0.5+i/HF_xy(1)/2,b,'Converting back to Time')
    end
    delete(b)
    
    
    %%%%% CONVOLUTION FILTERING %%%%%
if mode == 2
    
    if FsUS~=FsAE
        RefPulse       = resample(us,FsAE,FsUS);
    end
    RefPulse = RefPulse/(sum(abs(RefPulse)));
    
    %fprintf('Filtering 4D data\n')
    b = waitbar(0,'Filtering 4D data');
    %y = zeros(size(HF_xy,1),size(HF_xy,2),size(HF_zt(1))+length(RefPulse)-1,size(HF_zt(2)));
    
    Sz = HF_zt(1)+length(RefPulse)-1;
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            y{i,j} = zeros(Sz,HF_zt(2));
        end
    end
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            for k = 1:HF_zt(2)
                y{i,j}(:,k) = conv(HF{i,j}(:,k),RefPulse);
            end
        end
        waitbar(i/HF_xy(1),b,'Fast Time Filtering')
    end
    delete(b)
end
% M = 1.48*param.daq.HFdaq.pts/param.daq.HFdaq.fs_MHz;
% Mx = linspace(0,M,param.daq.HFdaq.pts);
% MxL = find(Mx < 50);
% MxH = find(Mx > 35);
% MxT = intersect(MxL,MxH);


% for i = 80:85
% test = real(y(:,1,MxT,i));
% test = squeeze(test);
% figure; imagesc(test')
% end


%vidwrite(real(y),MxT(1),MxT(end),1,size(y,1),20,120)


x = 2;
