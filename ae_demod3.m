function [DemodImage, sgnmtx, DemodImgdB] = ae_demod2(ModData,t,fc,maxref)
% function [DemodImage, sgnmtx, DemodImgdB] = ae_baseband(ModData,fs,fc,maxref)

narginchk(2,4);

[tLen, nCols] = size(ModData);

fs = 1/mean(diff(t));
if exist('fc','var') == 0 || isempty(fc) || isnan(fc)
    % determine the carrier frequency
    NFTL = max(2048,2^nextpow2(tLen));
    mod_s = double(mean(ModData,2));
    Y = fft(mod_s,NFTL);
    frq = linspace(0,1,NFTL/2+1)*fs/2;
    YAmp = abs(Y(1:length(frq)));
    
    % Set up fittype and options.
    ftObj = fittype( 'gauss1' );
    fitOpts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    fitOpts.Display = 'Off';
    fitOpts.Lower = [-Inf -Inf 0];
    fitOpts.StartPoint = [1 0.673828125 0.844670530491293];
    
    % Fit model to data.
    [fitres, gof] = fit( frq(:), YAmp(:), ftObj, fitOpts );
    fc = fitres.b1;
    
end

% 
% DemodImage = DemodFilt(1:tLen,:);
DemodImage = (-1) * amdemod(double(ModData),fc,fs);

% y = smooth(smooth(real(baseline)));
% end

sgnmtx = (-1)*sign(DemodImage);
% ModDataEnv = abs(hilbert(ModData));
% DemodImage = ModDataEnv.* sgnmtx;
ModDataEnv = abs(DemodImage);

if(nargout > 2)
    if(exist('maxref','var'))
        maxval = maxref;
    else
        maxval = max(ModDataEnv(:));
    end
    DemodImgdB = -20*log10(ModDataEnv/maxval + 1e-12).*sgnmtx;
end



