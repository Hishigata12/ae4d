% HF should be your 4D AE signal
% ave and in should contain 3 values, the first value should be 0 or 1 to
% tell the function whether to try that type of filter, the 2nd should be
% the x direction value, the third should be z direction value

function X = filts2D(HF,ave,in)
if ave(1) == 1
    Sx = zeros(size(HF));
    H = ones(ave(2),ave(3))/(ave(2)*ave(3));
    for i = 1:size(HF,4)
        for j = 1:size(HF,2)
            Sx(:,j,:,i) = filter2(H,squeeze(HF(:,j,:,i)));
        end
    end
    X = Sx;
else 
    X = HF;
end

if in(1) == 1
%     x = linspace(1,size(HF,1),size(HF,1));
%     z = linspace(0,size(HF,3),size(HF,3));
%     x2 = linspace(0,size(HF,1),size(HF,1)*in(2));
%     z2 = linspace(0,size(HF,3),size(HF,3)*in(3));
    
    [x, z] = meshgrid(1:size(HF,1),1:size(HF,3));
    [x2,z2] = meshgrid(1:(1/in(2)):size(HF,1),1:(1/in(3)):size(HF,3));
       for i = 1:size(HF,4)
        for j = 1:size(HF,2)
            Ix(:,j,:,i) = interp2(x,z,squeeze(HF(:,j,:,i))',x2,z2);
        end
       end
       X = permute(Ix,[3 2 1 4]);
end
            

