%x contains output filtered signal
%HF contains 4D AE data
%PE contains the US waveform
%mode signifies whether using match filter (1) or not (0)
%handles is the handles struct passed from GUI which contains tc info
%c is the cutoff frequencies (MHz) used for non-match-filter data

function y = w_ae_filt2(param,HF,PE,mode,handles,c)

%mode = 0;
%~~~~Create Filter~~~~%
tc = handles.tc.Value;
if tc
    tc_params = evalin('base','tc_params');
end
if mode == 1
    us = PE.TW.Wvfm1Wy;
    FsUS = PE.bScanParm.vsx_fs;
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
else
    Lae = param.daq.HFdaq.pts;
    FsAE = param.daq.HFdaq.fs_MHz;
    freqae = linspace(0,FsAE,Lae);
    max_fAE = length(freqae)/2;
end
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
    hsize = round(length(find(Hd == 1))/2);
    hd =ifft(Hd);
    if tc
        hwin = logspace(1,tc_params.freq_divisor,hsize)/10;
        if size(hwin,2) > size(hwin,1)
            hwin = hwin';
        end
        wq = hwin;
    else
        
        w1q = hamming(hsize);
        %     wq = circshift(w1q,round(hsize/2));
        %     h = padarray(wq,size(HF{1},1)-100,'post').*hd';
        %     H = fft(h,size(HF{1},1));
        wq = w1q;
    end
    H1 = padarray(wq,find(Hd == 1,1)-1,'pre');
    H2 = padarray(H1,length(Hd)-length(H1),'post');
    H2 = H2(1:end/2);
    H2 = [H2; flipud(H2)]; %cuz im a tricky motha fucka
    if length(H2) < length(Hd)
        zp = zeros(length(Hd) - length(H2),1);
        H2 = [zp; H2];
    end
    H = Hd.*H2';
    if mod(hsize,2) == 0
        m = 1;
    else 
        m = 0;
    end
    if tc
    H = imboxfilt(H,round(hsize/2)+m);
    end
    
       
    
    % elseif mode == 1
    %     if FsUS > FsAE
    %         RefPulse       = resample(us,FsAE,FsUS);
    %
    %     else
    %         RefPulse = resample(us,FsUS,FsAE);
    %     end
    %
    %     RefPulse = RefPulse/(sum(abs(RefPulse)));
    %     RefPulse = flipud(conj(RefPulse));
    %     H = fft(RefPulse);
    %     H = interp1(linspace(0,FsAE,length(RefPulse)),H',linspace(0,FsAE,Lae));
    %     a = find(freqae > 0.4,1);
    %     s = length(H);
    %     H2(1:(a-1)) = 0;
    %     H2(a:s-(a-1)) = 1;
    %     H2(s-(a-1):s) = 0;
    %     H = H.*H2;
    % end
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            % X{i,j} = HF{i,j}(1:end/2,:);
            X2{i,j} = zeros(HF_zt(1),HF_zt(2));
            y{i,j} = zeros(HF_zt);
        end
    end
    
    % fprintf('Filtering 4D data\n')
    %   b = waitbar(0,'Filtering 4D data');
    if size(H,2) > size(H,1)
        H = H';
    end
    if length(H) < size(HF{1,1},1)
        H = padarray(H,size(HF{1,1},1)-length(H),'pre');
    end
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            for k = 1:HF_zt(2)
                X2{i,j}(:,k) = fft(HF{i,j}(:,k)).*H;
                %X2(i,j,k,:) = filter(Hd,squeeze(X(i,j,k,:)));
            end
        end
        %  waitbar(i/HF_xy(1)/2,b,'Fast Time Filtering')
        multiWaitbar('Fast Time Filtering',i/HF_xy(1));
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
        %   waitbar(0.5+i/HF_xy(1)/2,b,'Converting back to Time')
        multiWaitbar('Converting to time domain',i/HF_xy(1));
    end
    % delete(b)
  x =2;  
    
    %%%%% CONVOLUTION FILTERING %%%%%
elseif mode == 1
    
    if FsUS~=FsAE
        RefPulse       = resample(us,FsAE,FsUS);
    end
    RefPulse = RefPulse/(sum(abs(RefPulse)));
    RefPulse = flipud(conj(RefPulse));
    %     hwin = hamming(length(RefPulse));
    %     RefPulse = hwin.*RefPulse;
    if tc
        hwin = logspace(tc_params.freq_divisor,1,length(RefPulse))/10;
        %hwin = linspace(0.3,1,length(RefPulse))';
        if size(hwin,2) > size(hwin,1)
            hwin = hwin';
        end
        RefPulse = RefPulse.*hwin;
    end
    
    %fprintf('Filtering 4D data\n')
    % b = waitbar(0,'Filtering 4D data');
    %y = zeros(size(HF_xy,1),size(HF_xy,2),size(HF_zt(1))+length(RefPulse)-1,size(HF_zt(2)));
    
    Sz = HF_zt(1)+length(RefPulse)-1;
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            y2{i,j} = zeros(Sz,HF_zt(2));
            y{i,j} = zeros(HF_zt(1),HF_zt(2));
        end
    end
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            for k = 1:HF_zt(2)
                y2{i,j}(:,k) = conv(HF{i,j}(:,k),RefPulse);
            end
        end
        %   waitbar(i/HF_xy(1),b,'Fast Time Filtering')
        multiWaitbar('Fast Time Filtering',i/HF_xy(1));
    end
    y = y2;
%     for i = 1:HF_xy(1)
%         for j = 1:HF_xy(2)
%             for k = 1:HF_zt(2)
%                 y{i,j}(:,k) = interp1(linspace(0,HF_zt(1),size(y2{1},1)),y2{i,j}(:,k),linspace(0,HF_zt(1),HF_zt(1)));
%             end
%         end
%         % waitbar(i/HF_xy(1),b,'Compressing Depth Axis');
%         multiWaitbar('Compressing Depth Axis',i/HF_xy(1));
%     end
    x=1;
    
    
    %  delete(b)
end



