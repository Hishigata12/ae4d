function wfilt = wieners(X,f,mu,b)
% x is the k-space data of the noise profile
%f is the frequency axis for the fft data

%b is the width of the expected signal (.2 is good for 1MHz)
%mu is the center frequency



% f_start = find(f>b(1),1);
% f_stop = find(f>b(2),1);
% f_mid = round((f_stop+f_start)/2);
% 
% s = linspace(f_start-f_start+1,f_mid-f_start,(f_mid-f_start))./(f_mid-f_start); 
% s2 = fliplr(s);
% s3 = [s s2(2:end)];

s = normpdf(mu,f,b); %generates gaussian waveform
%plot(f,s)
flen = length(f); %gets length of frequency vector
N_data = abs(X); %preps noise for normalization
s_norm = s./max(s); %normalizes gaus signal 
N_data_norm = N_data./max(N_data(round(flen/20):flen)); %normalizes noise
N_data_loc = find(max(N_data)); %gets location of max signal (whether positive or negative)
for i = 1:flen
    wfilt(i) =  s_norm(i)/(s_norm(i)+N_data_norm(i));
end
end
