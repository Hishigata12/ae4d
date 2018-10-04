%For R:
%   1 = xz
%   2 = yz
%   3 = xy
%   4 = zt

function [X,xaxis,yaxis] = plot_fft(param,img,R)
if R == 1
x_s = param.velmex.XDist/param.velmex.XNStep; %sampling rate for x axis
y_s = param.daq.HFdaq.fs_MHz; %sampling rate for y axis
dim = size(img); %dimensions of image
xaxis = linspace(0,x_s/2,round(dim(2)/2));
yaxis = linspace(0,y_s/2,round(dim(1)/2));
X = abs(fft(img));
elseif R == 4
x_s = param.daq.HFdaq.pulseRepRate_Hz; %sampling rate for x axis
y_s = param.daq.HFdaq.fs_MHz; %sampling rate for y axis
dim = size(img); %dimensions of image
xaxis = linspace(0,x_s/2,round(dim(2)/2));
yaxis = linspace(0,y_s/2,round(dim(1)/2));
X = abs(fft2(img));
end