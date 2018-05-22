% Outputs filter spatial domain projection data g(phi,s)
% Takes projection data p structured p(sum at s, angle number), absolute
% angle is not important

%Currently uses a symmetrical ramp filter to increase high frequencies, but
%other filters are also possible and probably better

%Version 1.0

%Author Chet Preston 2018

function q = cRadFilt(p)

n = size(p,1);
m = floor(n/2);

H(1:m) = linspace(0,1,m);
if mod(n,2) == 1
H(m+1:n) = linspace(1,0,m+1);
else
    H(m+1:n) = linspace(1,0,m);
end
    
H = H';

p = fft(p);

p = bsxfun(@times, p, H);

q = ifft(p,'symmetric');
%figure; imagesc(q);
assignin('base','Filter',H);
assignin('base','Filtered_Projections',q);
