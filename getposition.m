function B_mode_disp(x,t)
%x is 3dimensional matrix containing data
%t is a number or vector containing the time to view the B mode image

if length(t) == 1
t1 = x(:,:,t)
t2 = squeeze(t1)
imagesc(t2)
else
    
