function [xdemod3, X_bb1, t, f_ax, X_bb2, X_bb3] = baseband(X,fc,fs,wc1,wc2)
%X is frequency domain signal (after filtering)
%fc is center frequency
%fs is sampling frequency
%x_bb is basebanded time signal
%X_bb is basebanded freq signal
%t = linspace(0,size(X,1)/fs,size(X,1));

%t = (1:length(X)).'/fs;
%x_bb = X*sqrt(2).*exp(-1i*2*pi*fc*t);% Trying complex demod
%x_bb = X./cos(2*pi*fc*t); t cosine demod
N = size(X,1);
%half = floor(length(x_bb)/2);
NFFT = 2^nextpow2(N);
f_ax = fs/2*linspace(0,1,NFFT/2);
%f2 = fs/2*linspace(-1,1,NFFT); %Don't need since only need signal from
%0:pi
%X_bb = fft(x_bb,NFFT)/N;
%X_bb_shift = fftshift(x_bb,NFFT)/N;
%% For k > BB > t
%xmod = cos(2*pi*t*fc);
% xmod = sqrt(2).*exp(1i*2*pi*fc*t);
% xdemod2 = X.*xmod';
% Xmod = fft(xmod);
% Xfft = fft(xdemod2);
% for j = 1:size(X,2)
% X_bb(:,j) = conv(Xfft(:,j),Xmod);
% end
% X_bb1 = X_bb(1:length(f),:);
% X_bb2 = X_bb(length(f)-1:end,:);
% X_bb3 = [X_bb2(1:end-1,:);X_bb1(1:end-1,:)];
% t_ax = linspace(0,length(X)/fs,length(X));
% x_bb = ifft(X_bb1);
% f_ax = f;
if length(size(X)) > 2
    mode = 1;
else
    mode = 2;
end

%Gets analytic signal from Hilbert Transform xa = x+jxh
m = 1j;
dims = size(X);
%t = (1:dims(3)) /fs;
if length(dims) < 4
    dims(4) = 1;
end

% h = 1./(pi.*t);
% b = waitbar(0,'Basebanding');
if mode == 1
    for i = 1:dims(1)
        for j = 1:dims(2)
            for k = 1:dims(4)
                %             H(i,j,:,k) = conv(squeeze(X(i,j,:,k)),h');
                x(i,j,:,k) = fft(X(i,j,:,k));
            end
        end
%         waitbar(i/dims(1)/4,b,'Computing FFT')
multiWaitbar('Computing FFT',i/dims(1));
    end
elseif mode  == 2
    for i = 1:dims(1)
        x(i,:) = fft(X(i,:));
        multiWaitbar('Computing FFT',i/dims(1));
    end
end
% S = -m*sign(x);
% xh = x+S;
% x1 = x + xh;

% xhilb = hilbert(X);
% xa1 = X + m*xhilb;
% for i = 1:dims(1)
%     for j = 1:dims(2)
%         for k = 1:dims(4)
%             xdemod(i,j,:,k) = squeeze(xa1(i,j,:,k)).*exp(-m*2*pi*t'*fc);
%         end
%     end
% end



if mode == 1
    x2 = 2*x(:,:,1:end/2,:);
    f_axis = linspace(0,fs/2,size(x2,3));
    wlow = find(f_axis > wc1,1)-1;
    whigh = find(f_axis >= wc2,1);
    hamlen = whigh-wlow;
    Hwin = hamming(hamlen);
    Hwin = padarray(Hwin,wlow-1,0,'pre');
    Hwin = padarray(Hwin,size(x2,3)-whigh+1,0,'post');
    f_ax = fs/2*linspace(0,1,dims(3)/2);
elseif mode == 2
    x2 = 2*x(:,1:end/2);
    f_axis = linspace(0,fs/2,size(x2,2));
    wlow = find(f_axis > wc1,1)-1;
    whigh = find(f_axis >= wc2,1);
    hamlen = whigh-wlow;
    Hwin = hamming(hamlen);
    Hwin = padarray(Hwin,wlow-1,0,'pre');
    Hwin = padarray(Hwin,size(x2,2)-whigh+1,0,'post');
    f_ax = fs/2*linspace(0,1,dims(2)/2);
end

%x2(:,:,1:2,:) = 0;



% for i = 1:dims(1)
%     for j = 1:dims(2)
%         for k = 1:dims(4)
%             x3(i,j,:,k) = fft(H(i,j,:,k));
%         end
%     end
% end


% a = squeeze(x1(15,1,:,21));
% a2 = squeeze(x2(15,1,:,21));
% a3 = squeeze(x3(15,1,:,21));

% figure; plot(abs(a))
%figure; plot(f_ax,abs(a2));
% figure; plot(abs(a3))


dims = size(x2);
if length(dims) < 4
    dims(4) = 1;
end

pk = zeros(size(x2,1),size(x2,2));




if mode == 1
    for i = 1:dims(1)
        for j = 1:dims(2)
            for k = 1:dims(4)
                x2(i,j,:,k) = Hwin.*squeeze(x2(i,j,:,k));
                X2(i,j,:,k) = ifft(x2(i,j,:,k));
                if max(x2(i,j,:,k)) ~= 0
                    pk(i,j,k) = f_ax(find(abs(x2(i,j,:,k)) == max(abs(x2(i,j,:,k)))));
                else
                    pk(i,j,k) = 0;
                end
            end
            
        end
%         waitbar(.25 + i/dims(1)/4,b,'Filtering')
multiWaitbar('Filtering',i/dims(1));
    end
elseif mode == 2
    for i = 1:dims(1)
        x2(i,:) = Hwin'.*squeeze(x2(i,:));
        X2(i,:) = ifft(x2(i,:));
        if max(x2(i,:)) ~=0
            pk(i) = f_ax(find(abs(x2(i,:)) == max(abs(x2(i,:)))));
        else
            pk(i) = 0;
        end
%         waitbar(.25 + i/dims(1)/4,b,'Filtering')
multiWaitbar('Computing FFT',i/dims(1));
    end
    
end


if fc == 0
    fc = pk;
else
    fc2 = ones(size(pk));
    fc = fc2.*fc;
end

if mode == 1
    t = (1:dims(3)) /fs;
    for i = 1:dims(1)
        for j = 1:dims(2)
            for k = 1:dims(4)
                xdemod2(i,j,:,k) = squeeze(X2(i,j,:,k)).*exp(-m*2*pi*fc(i,j,k)*t');
                %xdemod2(i,j,:,k) = squeeze(X(i,j,:,k)).*exp(-m*2*pi*fc*t');
                xdemod3(i,j,:,k) = interp1(linspace(0,1,dims(3)),squeeze(xdemod2(i,j,:,k)),linspace(0,1,size(X,3)));
            end
        end
%         waitbar(.5 + i/dims(1)/2,b,'Demodulating')
multiWaitbar('Demodulating',i/dims(1));
    end
elseif mode == 2
    t = (1:dims(2))/fs;
    for i = 1:dims(1)
        xdemod2(i,:) = squeeze(X2(i,:))'.*exp(-m*2*pi*fc(i)*t');
        xdemod3(i,:) = interp1(linspace(0,1,dims(2)),xdemod2(i,:),linspace(0,1,size(X,2)));
%         waitbar(.5 + i/dims(1)/2,b,'Demodulating')
        multiWaitbar('Demodulating',i/dims(1));
    end
end

% for i = 1:dims(1)
%     for j = 1:dims(2)
%         for k = 1:dims(4)
%            % xdemod3(i,j,:,k) = interp1(linspace(0,1,dims(3)),squeeze(xdemod2(i,j,:,k)),linspace(0,1,size(X,3)));
%         end
%     end
%     waitbar(.75 + i/dims(1)/4,b,'Interpolating')
% end

% delete(b)
%figure; imagesc(real(squeeze(xdemod3(:,1,:,21))));
% figure; imagesc(real(squeeze(xdemod(:,1,:,21))));
%
% for i = 1:dims(1)
%     for j = 1:dims(2)
%         for k = 1:dims(4)
%             Xdemod2(i,j,:,k) = fft(squeeze(xdemod3(i,j,:,k)));
%             Xdemod(i,j,:,k) = fft(squeeze(xdemod(i,j,:,k)));
%         end
%     end
% end

%figure; plot(abs(squeeze(Xdemod2(15,1,:,21))));
% figure; plot(abs(squeeze(Xdemod(15,1,:,21))));

%% For t ==> BB
% dims = size(X);
% t = (1:dims(3)) /fs;
% for i = 1:dims(1)
%     for j = 1:dims(2)
%         if length(dims) == 4
%             for k = 1:dims(4)
%                 xdemod2(i,j,:,k) = squeeze(X(i,j,:,k)).*sqrt(2).*exp(-1i*2*pi*fc*t)';
%             end
%         else
%             xdemod2(i,j,:) = squeeze(X(i,j,:)).*sqrt(2).*exp(-1i*2*pi*fc*t)';
%         end
%     end
% end

%xdemod2 = X.*xmod';
% Xmod = fft(xmod);
% X_bb1 = fft(xdemod2);

%% for BB > k> LPF > t using sinusoidal demod
% xmod = cos(2*pi*t*fc);
% xdemod = xmod'./X;
% Xf = fft(X,NFFT);
% Xdemod = fft(xdemod,NFFT);
% LPF_loc = find(f > 2,1);
% LPF(1:LPF_loc) =1;
% LPF(LPF_loc + 1: size(Xdemod,1)) = 0;
% Xlpf = Xdemod.*LPF';
% x_lpf = ifft(Xlpf);

% imagesc(m,Ax,real(X));
% center_axis(real(X),Ax,'HF',[45 65],'lin');
% figure
% imagesc(m,Ax,abs(xdemod2));
% center_axis(abs(xdemod2),Ax,'HF',[45 65],'lin');

% Design LP filter
if length(size(xdemod3)) == 2
    xdemod3 = permute(xdemod3,[1 3 2]);
end


hd = designfilt('lowpassfir','PassbandFrequency',0.15, ...
         'StopbandFrequency',0.2,'PassbandRipple',0.5, ...
         'StopbandAttenuation',30,'DesignMethod','kaiserwin');
     for i = 1:size(xdemod3,1)
         for j = 1:size(xdemod3,2)
             for k = 1:size(xdemod3,4)
                 xdemod(i,j,:,k) = interp1(linspace(0,1,size(xdemod3,3)),squeeze(xdemod3(i,j,:,k)),linspace(0,1,size(xdemod3,3)*2));
        xdemod4(i,j,:,k) = filtfilt(hd,squeeze(xdemod(i,j,:,k)));
        xdemod1(i,j,:,k) = interp1(linspace(0,1,size(xdemod3,3)),squeeze(xdemod3(i,j,:,k)),linspace(0,1,size(xdemod,3)/2));
        snips = round(size(xdemod1,3)*0.15);
        trunct = xdemod1(i,j,snips:(size(xdemod1,3)-snips),k);
        padded = padarray(squeeze(trunct),snips,'symmetric','both');
        xdemod5(i,j,:,k) = padded(1:size(xdemod1,3));
        multiWaitbar('LPF',i/size(xdemod3,1));
             end
         end
     end
     xdemod3 = xdemod5;

end