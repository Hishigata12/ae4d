% HF should be your 4D AE signal
% ave and in should contain 3 values, the first value should be 0 or 1 to
% tell the function whether to try that type of filter, the 2nd should be
% the x direction value, the third should be z direction value

function X = filts3D(HF,ave,in,med,p,param)
%dims = [param.velmex.XNStep param.velmex.YNStep param.daq.HFdaq.pts param.daq.HFdaq.NoBurstTriggers];
dims = size(HF);
param = param;
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
        multiWaitbar('Median Filtering',i/dims(4));
    end
    HF = Sx;
    if p(1) > 1
        x = HF;
        if mod(p(1),2) == 0
            errordlg('window must be odd')
        else
            for i = 1:size(x,1)
                for j = 1:size(x,2)
                    for k = 1:size(x,3)
                        x2(i,j,k,:) = medfilt1(squeeze(x(i,j,k,:)),p(1));
                    end
                end
                multiWaitbar('Median Filtering in Time',i/size(x,1));
            end
        end
        HF = x2;
        clear x2;
    end
end

if ave(1) == 1
    HF = real(HF);
    dims = size(HF);
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
        Sx(:,:,:,i) = imboxfilt3(HF(:,:,:,i),ave(2:end));
        multiWaitbar('Smoothing Signal',i/dims(4));
    end
    clear HF
    HF = Sx;
    if p(2) > 1
        x = HF;
        if mod(p(2),2) == 0
            errordlg('window must be odd')
        else
            a = 1;
            b = (1/p(2))*ones(1,p(2));
            for i = 1:size(x,1)
                for j = 1:size(x,2)
                    for k = 1:size(x,3)
                        x2(i,j,k,:) = filter(b,a,squeeze(x(i,j,k,:)));
                    end
                end
                multiWaitbar('Smoothing in Time',i/size(x,1));
            end
        end
        HF = x2;
        clear x2;
    end
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
        q = sum(P);
        
    end
    if q == 0
        [x, y, z] = meshgrid(1:dims(2),1:dims(1),1:dims(3));
        [x2, y2, z2] = meshgrid(1:(1/in(3)):dims(2),1:(1/in(2)):dims(1),1:(1/in(4)):dims(3));
        %Ix = zeros(size(x2));
        for i = 1:dims(4)
            if in(5) == 0
                Ix(:,:,:,i) = interp3(x,y,z,HF(:,:,:,i),x2,y2,z2);
            elseif in(5) == 1
                Ix(:,:,:,i) = squarify(HF(:,:,:,i),'m');
            end
            multiWaitbar('Interpolating',i/dims(4))
        end
        X = Ix;
        
    end
    if q > 0
        in = in([1 2 4 5]);
        [x, z] = meshgrid(1:dims(1),1:dims(3));
        [x2,z2] = meshgrid(1:(1/in(2)):dims(1),1:(1/in(3)):dims(3));
        for i = 1:dims(4)
            for j = 1:dims(2)
                if in(4) == 0
                    Ix(:,j,:,i) = interp2(x,z,squeeze(HF(:,j,:,i))',x2,z2);
                elseif in(4) == 1
                    Ix(:,j,:,i) = permute(squarify(HF(:,j,:,i),'m'),[1 3 2]);
                end
                multiWaitbar('Interpolating',i/dims(4));
            end
        end
        
        X = permute(Ix,[3 2 1 4]);
    end
    if p(3) < 1
        Q = 'Decimating';
    elseif p(3) > 1
        Q = 'Interpolating';
    end
    clear x2
    if p(3) ~= 1
        x = X;
        
        for i = 1:size(x,1)
            for j = 1:size(x,2)
                for k = 1:size(x,3)
                    x2(i,j,k,:) = interp1(linspace(0,1,size(x,4)),squeeze(x(i,j,k,:)),linspace(0,1,size(x,4)*p(3)));
                end
            end
            multiWaitbar([Q ' in Time'],i/size(x,1));
        end
        X = x2;
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

multiWaitbar('CLOSEALL');


