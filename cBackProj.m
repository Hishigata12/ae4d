%Author: Chet Preston
%2018
%Adapted from iradon.m
% IM2: sinogram with dimension 1 = line integrals and dimension 2 =
% angles
% A: range of theta in radians
% f: 1 = centered at min, 2 = centered at middle


function [BPI] = cBackProj(IM2,A,f,w)

n = size(IM2,1);
m = size(IM2,2);
Theta = linspace(0,A,m); %Pi will need to be variable based on range of angles used; possibly pi/4
if f == 2
    Theta = Theta-median(Theta); %Changes center location
end

for i = 1:size(IM2,2)
    G(:,i) = fft(IM2(:,i));
end

%%%%%%%%Filter the projections here%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BPI = zeros(w,n);

midpoint = ceil(n/2);


% Final reconstructed discrete image will follow:
% G(w,phi)*exp(-jwm(sin(phi))/M+jwn(cos(phi))/N)

x = 1:n;
y = 1:w;
[X,Y] = meshgrid(x,y);
xpr = X - midpoint -1;
ypr = Y - midpoint -1;

ct = cos(Theta);
st = sin(Theta);


for i = 1:m
    proj = IM2(:,i); % gives filtered line of projection data at given theta
    taxis = (1:size(IM2,1)) - midpoint;
    t = xpr.*st(i)+ypr.*ct(i); % Might need to flip ct and st. x and y might need changing.
    projContrib = interp1(taxis,proj,t(:),'linear');
    BPI = BPI + reshape(projContrib,w,n);
end

BPI = BPI*pi/(2*length(Theta));
    
