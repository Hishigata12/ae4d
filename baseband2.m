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

%Gets analytic signal from Hilbert Transform xa = x+jxh
m = 1j;
dims = size(X);
t = (1:dims(3)) /fs;
if length(dims) < 4
    dims(4) = 1;
end

% h = 1./(pi.*t);
b = waitbar(0,'Basebanding');
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(4)     
%             H(i,j,:,k) = conv(squeeze(X(i,j,:,k)),h');
            x(i,j,:,k) = fft(X(i,j,:,k));
        end
    end
    waitbar(i/dims(1)/4,b,'Computing FFT')
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




x2 = 2*x(:,:,1:end/2,:);
f_axis = linspace(0,fs,size(x2,3));
wlow = find(f_axis > wc1,1)-1;
whigh = find(f_axis >= wc2,1);
hamlen = whigh-wlow;
Hwin = hamming(hamlen);
Hwin = padarray(Hwin,wlow-1,0,'pre');
Hwin = padarray(Hwin,size(x2,3)-whigh+1,0,'post');

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
f_ax = fs/2*linspace(0,1,dims(3)/2);
% figure; plot(abs(a))
%figure; plot(f_ax,abs(a2));
% figure; plot(abs(a3))


dims = size(x2);
if length(dims) < 4
    dims(4) = 1;
end


for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(4)
            x2(i,j,:,k) = Hwin.*squeeze(x2(i,j,:,k));
            X2(i,j,:,k) = ifft(x2(i,j,:,k));
        end
    end
     waitbar(.25 + i/dims(1)/4,b,'Filtering')
end
t = (1:dims(3)) /fs;
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(4)
            xdemod2(i,j,:,k) = squeeze(X2(i,j,:,k)).*exp(-m*2*pi*fc*t');
             %xdemod2(i,j,:,k) = squeeze(X(i,j,:,k)).*exp(-m*2*pi*fc*t');
        end
    end
       waitbar(.5 + i/dims(1)/4,b,'Demodulating')
end
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(4)
            xdemod3(i,j,:,k) = interp1(linspace(0,1,dims(3)),squeeze(xdemod2(i,j,:,k)),linspace(0,1,size(X,3)));
        end
    end
    waitbar(.75 + i/dims(1)/4,b,'Interpolating')
end

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

%Design LP filter

end