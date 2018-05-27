function [xdemod2, X_bb1, t, f_ax, X_bb2, X_bb3] = baseband(X,fc,fs)
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

%% For t > BB
dims = size(X);
t = (1:dims(3)) /fs;
for i = 1:dims(1)
    for j = 1:dims(2)
        if length(dims) == 4
            for k = 1:dims(4)
                xdemod2(i,j,:,k) = squeeze(X(i,j,:,k)).*sqrt(2).*exp(-1i*2*pi*fc*t)';
            end
        else
            xdemod2(i,j,:) = squeeze(X(i,j,:)).*sqrt(2).*exp(-1i*2*pi*fc*t)';
        end
    end
end

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