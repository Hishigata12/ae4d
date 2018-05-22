function IM2 = cRadon(IM1,n)

% pad image with zeros if rectangular 

if size(IM1,3) > 1
    IM1 = rgb2gray(IM1);
end

if isinteger(IM1)
    IM1 = im2double(IM1); % Converts to double if originally integer
end

% [Xq,Yq] = meshgrid(0:1/d:size(IM1,1),0:1/d:size(IM1,2));
% Xq = Xq(1,:);
% Yq = Yq(:,1);


% if d > 1 
%     IM1 = interp2(IM1,Xq,Yq);
% elseif d < 1
%     IM1 = interp2(IM1,Xq,Yq);
% else
%     IM1 = IM1;
% end

[lIM, wIM] = size(IM1); %Gets size of image to determine needed padding

iDiag = sqrt(lIM^2 + wIM^2);
wPad = ceil(iDiag - wIM) +2; %Creates padded dimensions
lPad = ceil(iDiag - lIM) +2; 
padIM = zeros(lIM+lPad,wIM+wPad); %Creates larger padded array for image to rotate inside
padIM(ceil(lPad/2):(ceil(lPad/2)+lIM-1),ceil(wPad/2):ceil(wPad/2)+wIM-1) = IM1; %Centers image in padded array

%Loops through angles of images to convert 2D image f(x,y) into 1D sinogram of
%theta and length g(phi,s)

Theta = linspace(0,180,n); %creates array of different angles
IM2 = zeros(size(padIM,2),n); %creates array g(phi,s) which will contain the summed projections
for i = 1:n
    tmpIM = imrotate(padIM, Theta(i), 'bilinear', 'crop'); %performs pseudo angle (phi) shifts by rotating image
    IM2(:,i) = sum(tmpIM); %takes line integral (sum) at each point (s) along the rotated image
                          % Note that for this simplified discrete setup the lengths
                          % of s, x, and y are the same.
end
figure; imshow(IM2);
assignin('base','Projections',IM2);


