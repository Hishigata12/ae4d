function [x, X_filt2] = guassfilt(X,f,w)

L = length(X)/2+1;
s = normpdf(0,f,w);
X_filt = s.*X';
X_filt2 = X_filt';
x = ifft(X_filt2,length(f),1)*length(f);



end