% HF should be your 4D AE signal
% ave and in should contain 3 values, the first value should be 0 or 1 to
% tell the function whether to try that type of filter, the 2nd should be
% the x direction value, the third should be z direction value

function X = filts3D(HF,ave,in,med,param)
 b = waitbar(0);
 %dims = [param.velmex.XNStep param.velmex.YNStep param.daq.HFdaq.pts param.daq.HFdaq.NoBurstTriggers];
 dims = size(HF);
if length(dims) ==3
    dims(4) = 1;
end

if med(1) == 1
 HF = real(HF);
    dims = size(HF);
    for i = 2:4
        if mod(med(i),2) ~= 1
            errordlg('Need window size for each dimension to be odd')
            return
        end
    end
    Sx = zeros(dims);
    if length(dims) < 4
        dims(4) = 1;
    end
    % H = ones(ave(2),ave(3))/(ave(2)*ave(3));
  
    for i = 1:dims(4)
        Sx(:,:,:,i) = medfilt3(HF(:,:,:,i),med(2:end));
        waitbar(i/dims(4),b,'Median Filtering');
    end
    HF = Sx;
end


if in(1) == 1
    %     x = linspace(1,size(HF,1),size(HF,1));
    %     z = linspace(0,size(HF,3),size(HF,3));
    %     x2 = linspace(0,size(HF,1),size(HF,1)*in(2));
    %     z2 = linspace(0,size(HF,3),size(HF,3)*in(3));
    for i = 1:3
        if size(HF,i) == 1
            P(i) = 1;
        else
            P(i) = 0;
        end
        p = sum(P);
    end
    if p == 0 
            [x, y, z] = meshgrid(1:dims(2),1:dims(1),1:dims(3));
            [x2, y2, z2] = meshgrid(1:(1/in(3)):dims(2),1:(1/in(2)):dims(1),1:(1/in(4)):dims(3));
            for i = 1:dims(4)
                Ix(:,:,:,i) = interp3(x,y,z,HF(:,:,:,i),x2,y2,z2);
                waitbar(i/dims(4),b,'Interpolating')
            end
            X = Ix;
            
    end
    if p > 0
        in = in([1 2 4]);
        [x, z] = meshgrid(1:dims(1),1:dims(3));
    [x2,z2] = meshgrid(1:(1/in(2)):dims(1),1:(1/in(3)):dims(3));
       for i = 1:dims(4)
        for j = 1:dims(2)
            Ix(:,j,:,i) = interp2(x,z,squeeze(HF(:,j,:,i))',x2,z2);
             waitbar(i/dims(4),b,'Interpolating')
        end
       end
       X = permute(Ix,[3 2 1 4]);   
    end
%       b1 = round(dims(3)/10);
%     b2 = round(dims(3)*0.9);
%     X = real(20*log10(real(X)./max(max(max(real(X(:,:,b1:b2,:)))))));
else
    X = HF;
%     b1 = round(dims(3)/10);
%     b2 = round(dims(3)*0.9);
%     X = real(20*log10(real(X)./max(max(max(real(X(:,:,b1:b2,:)))))));
end
if ave(1) == 1
    X = real(X);
    dims = size(X);
    for i = 2:4
        if mod(ave(i),2) ~= 1
            errordlg('Need window size for each dimension to be odd')
            return
        end
    end
    Sx = zeros(dims);
    if length(dims) < 4
        dims(4) = 1;
    end
    % H = ones(ave(2),ave(3))/(ave(2)*ave(3));
  
    for i = 1:dims(4)
        Sx(:,:,:,i) = imboxfilt3(X(:,:,:,i),ave(2:end));
        waitbar(i/dims(4),b,'Smoothing Signal');
    end
    X = Sx;
    delete(b)
else
   % X = HF;
    delete(b)
end
        
        
