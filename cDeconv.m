%This function performs deblurring and denoising with Wiener filter
% b is psf 
% d is the depth range of signal
% v is time point of noise


function X = cDeconv(HF,b,d,v)
Y = squeeze(HF(:,1,:,v));
 n_var = var(Y(:));
b = normpdf(linspace(-2,2,100),0,.2);
for i = 1:size(HF,2)
    for j = 1:size(HF,4)
        Y = squeeze(HF(:,i,d,j));
        s_var = var(Y(:));
        X(:,i,:,j) = deconvwnr(squeeze(HF(:,i,:,j)),b,n_var/s_var);
    end
end