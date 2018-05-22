function R = cFlashRecon(b,k) % b is png file, k is number of angles


c = cRadon(b,k); % Performs RT

% cRadon will convert image into a squarish matrix and we'll need to
% convert back to native lengths later
crops = length(c) - size(b); % gets added vector sizes from cRadon
wtail = ceil(crops(1)/2); % Can add if statements to make correct adjustments to length for even length arrays
ltail = ceil(crops(2)/2);

q1 = cRadFilt(c);  % Change filter around for performance     

R = cBackProj(q1,pi,1); % Performs back projection operation to original image

R = fliplr(imrotate(R,-90)); % rotates image back to native direction

R = R(wtail:size(R,1)-wtail,ltail:size(R,2)-ltail); %crops image to original dimensions

end