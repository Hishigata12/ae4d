%x contains output filtered signal
%HF contains 4D AE data
%LF contains the injected current waveform
%n is the point along fast time point that contains the noise spectrum in the
%current direction

function [y lf] = w_slow_filt2(param,HF,LF,mode,c)
%~~~~Create Filter~~~~%

%mode = 0;

if length(size(HF{1}))>2
    for i = size(HF,1)
        for j = size(HF,2)
            HF2{i,j} = squeeze(HF{i,j}(:,:,2)); % Gets only a single channel for AE
        end
    end
    clear HF
    HF = HF2;
    clear HF2
end

% if size(LF,2) > 1
%     LF = LF(:,1); %Not sure why I had this in...
% end


Fs = param.daq.LFdaq.fs_Hz;
L = param.daq.LFdaq.pts;
freq = linspace(0,Fs,L);
Fs_us = param.daq.HFdaq.pulseRepRate_Hz;


L_us = param.daq.HFdaq.NoBurstTriggers;
f_us = linspace(0,Fs_us,L_us);

max_f = length(freq)/2;
max_us = L_us/2;



% d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',100/Fs_us*2,150/Fs_us*2,250/Fs_us*2,300/Fs_us*2,60,1,60);
% Hd = design(d,'equiripple');

%~~~~~Filter in slow time~~~~%
fspec = '* %d percent complete \n*';
HF_xy = size(HF);
HF_zt = size(HF{1});


%%%%% FREQUENCY FILTERING %%%%%%%%%%%%

if mode == 0
    
    Fc1 = find(f_us > c(1),1);%/Fs_us*2;
    Fc2 = find(f_us > c(2),1);%Fs_us*2;
    
    for i=1:L_us
        if (i >= Fc1 && i <= Fc2) || (i >= L_us-Fc2 && i <= L_us-Fc1)
            S_Filter(i) = 1;
        else
            S_Filter(i) = 0;
        end
    end

    s_filter = ifft(S_Filter);
    W = hamming(L_us);
    W2 = circshift(W,max_us);
    h = W2.*s_filter';
    H = fft(h);
    Hs = H(1:end/2);
    
% elseif mode == 1
%     if Fs > Fs_us
%         RefPulse       = resample(LF,Fs_us,Fs);
%    
%     elseif Fs < Fs_us
%         RefPulse = resample(LF,Fs,Fs_us);
%     end
%     
%     RefPulse = RefPulse/(sum(abs(RefPulse)));
%     H = fft(RefPulse);
%     a = find(f_us > 50,1);
%     s = length(RefPulse);
%     H2(1:(a-1)) = 0;
%     H2(a:s-(a-1)) = 1;
%     H2(s-(a-1):s) = 0;
%     H = H.*H2';
% end
%     
        
        
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            % X{i,j} = HF{i,j}(:,1:end/2);
            X2{i,j} = zeros(HF_zt(1),HF_zt(2));
            y{i,j} = zeros(HF_zt);
        end
    end
    
 %   fprintf('Filtering 4D data\n')
    %b = waitbar(0,'Filtering 4D data');
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            for k = 1:HF_zt(1)
                y{i,j}(k,:) = real(ifft(fft(HF{i,j}(k,:)).*H',HF_zt(2)));
                %X2(i,j,k,:) = filter(Hd,squeeze(X(i,j,k,:)));
                
            end
        end
        %waitbar(i/HF_xy(1)/2,b,'Slow Time Filtering')
        multiWaitbar('Slow Time Filtering',i/HF_xy(1));
    end
    
%     
    LF2 = fft(LF);
    LFHam = hamming(length(LF));
    lf_axis = linspace(0,Fs,L);
    Fc1 = find(lf_axis >= c(1),1);
    Fc2 = find(lf_axis >= c(2),1);
    for i=1:L
        if (i >= Fc1 && i <= Fc2) || (i >= L-Fc2 && i <= L-Fc1)
            LF_Filter(i) = 1;
        else
            LF_Filter(i) = 0;
        end
    end
    LFtime = ifft(LF_Filter);
   % LFH = interp1(linspace(1,10,L_us),H,linspace(1,10,L));
   LFt2 = circshift(LFtime,round(L/2));
   LFH1 = LFt2.*LFHam'; 
   LFH2 = fft(LFH1);
    for i = 1:size(LF,2)
        LF3(:,i) = LF2(:,i).*LFH2';
        lf(:,i) = real(ifft(LF3(:,i)));
    end
    lf = circshift(lf,round(L/2));
%     LF3 = LF2.*LFH';
%     lf = real(ifft(LF3));
% %    fprintf('Converting to space domain\n');
%     
%     %Xf = flip(X2,4);
%     %Xs = cat(4,X2,Xf);
%     
%     %y = zeros(size(Xs));
%     
%     for i = 1:HF_xy(1)
%         for j = 1:HF_xy(2)
%             for k = 1:HF_zt(1)
%                 y{i,j}(k,:) = ifft(X2{i,j}(k,:),HF_zt(2));
%             end
%         end
%       %  waitbar(0.5+i/HF_xy(1)/2,b,'Slow Time Filtering')
%         multiWaitbar('Converting back to time domain',i/HF_xy(1));
%     end
%  %   delete(b)
    
elseif mode == 1
    
    %%%%%CONVOLUTION FILTERING %%%%%%%%%%%%%
    
    if Fs~=Fs_us
        RefPulse       = resample(LF(:,1),Fs_us,Fs);
    end
    
    RefPulse = RefPulse/(sum(abs(RefPulse)));
  %  RefPulse = flipud(conj(RefPulse));
    
%    fprintf('Filtering 4D data\n')
   % b = waitbar(0,'Filtering 4D data');
    %y = zeros(size(HF_xy,1),size(HF_xy,2),size(HF_zt(1))+length(RefPulse)-1,size(HF_zt(2)));
    
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            y{i,j} = zeros(HF_zt(1),HF_zt(2)+length(RefPulse)-1);   
            Q{i,j} = zeros(HF_zt(1),HF_zt(2)+length(RefPulse)-1);   
        end
    end
    
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            for k = 1:HF_zt(1)
                y{i,j}(k,:) = conv(HF{i,j}(k,:),RefPulse);
            end
        end
       % waitbar(i/HF_xy(1),b,'Slow Time Filtering')
       multiWaitbar('Slow Time Filtering',i/HF_xy(1));
    end
    %
    %Removes DC components
    Hfilt = ones(1,size(y{1,1},2));
    Hfilt(1:2) = 0;
    Hfilt(end-1:end) = 0;
    for i = 1:HF_xy(1)
        for j = 1:HF_xy(2)
            for k = 1:size(y{1,1},1)
                Q{i,j}(k,:) = fft(y{i,j}(k,:));
                Q{i,j}(k,:) = Q{i,j}(k,:).*Hfilt;
                y{i,j}(k,:) = ifft(Q{i,j}(k,:));
            end
        end
         multiWaitbar('Removing DC',i/HF_xy(1));
    end
  %  delete(b)
  lf = LF;
end

% for i = 925
%     figure;plot(real(y{13,1}(i,:)))
% end



