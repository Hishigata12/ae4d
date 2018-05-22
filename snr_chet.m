function db = snr_chet(X, xsig, xnoise)

%X = full AE signal
%xsig = bounds of AE signal
%xnoise = bounds of equiproportional noise signal

s = X(xsig(1):xsig(2));
n = X(xnoise(1):xnoise(2));
srms = rms(s);
snoise = rms(n);
db = 20*log10(srms/snoise);