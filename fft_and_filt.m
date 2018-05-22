function [xfilt, Xf2, f, filt_win] = fft_and_filt(x, Fs, filt_type, b,c, w)

global m_burst_pts
global X
%perform FFT on data
%x is the time domain data
%Fs = param.daq.LFdaq.fs_Hz;
%filt_type denotes filter: can be 'Hann' 'Wiener' or 'none';
%b is window for hann or burst number for wiener
%c is center frequency for wiener filter
%w is variance of wiener filter


% if exist('b','var')
% if length(b) ~= 4 && filt_type(1:4) == 'Hann'
%     errordlg('b must be ascending vector of length 4','Bad filter window')
% elseif length(b) ~= 1 && filt_type(1:4) == 'Wien'
%     errordlg('b must be single number','Bad filter window')
% end
% end

if filt_type(1:4) == 'Wien'
    xn = x(:,b);
    x2 = x;
    clear x;
    x = xn;
    if size(x2,1) > size(x2,2)
    x_red = x(600:end);  %gets rid of the big bang **assumes 20MHz Fs** Starts at 20mm
    x_pad = padarray(x_red,599,0,'pre');
    N = size(x_pad,1);
    else
        x_red = x(ceil(size(x,1))/1.6:end);
        x_pad = padarray(x_red,floor(size(x,1))/2,'post');
        N = size(x_pad,1);
    end
NFFT = 2^nextpow2(N);
X = fft(x_pad,NFFT)/N;
f = Fs/2*linspace(0,1,NFFT/2+1);
flen = length(f);
jump = f(5)-f(4);
X2 = fft(x2,NFFT)/N;
else   
N = size(x,1);
NFFT = 2^nextpow2(N);
X = fft(x,NFFT)/N;
f = Fs/2*linspace(0,1,NFFT/2+1);
flen = length(f);
jump = f(5)-f(4);
end



%Generate Hann window filter
if filt_type(1:4) == 'Hann'
for i = 1:4
    bloc(i) = find(f>b(i),1);
    
end
filt_win(1:bloc(1)) = 0;
slope1 = bloc(2)-bloc(1);%round((b(1)+lf_jump:b(2))/lf_jump);
slope2 = bloc(4)-bloc(3);%round((b(3)+lf_jump:b(4))/lf_jump);
filt_win(bloc(1)+1:bloc(2)) = (1./(slope1./(1:slope1))); %(1-sqrt(2:slope1)))/(1/(1-sqrt(slope1)));
filt_win(bloc(2)+1:bloc(3)) = 1;
filt_win(bloc(3)+1:bloc(4)) = (1-(1:slope2).^2)/slope2^2+1;
filt_win(bloc(4)+1:length(f)) = 0;
elseif filt_type(1:4) == 'none'
elseif filt_type(1:4) == 'Wien'
    wfilt = wieners(X,f,c,w);
else
    errordlg('filter type must be either Hann, Weiner or none','unknown or incorrect filter') 
end

%Filter k space
if filt_type(1:4) == 'Hann'
    % for i = 1:size(X,2)
    Xf = X(1:flen,:).*filt_win';
    %    end
elseif filt_type(1:4) == 'Wien'
    
    Xf = X2(1:flen,:).*wfilt';
else
    Xf = X;
end



%Inverse FFT
if size(Xf,1) < size(Xf,2)
    
    xfilt = ifft(Xf,NFFT)*N;
    xfilt = xfilt';
   % xfilt = flipud(xfilt);
    Xf2 = Xf';
else
    xfilt = ifft(Xf,NFFT)*N;
    Xf2 = Xf;
end

end