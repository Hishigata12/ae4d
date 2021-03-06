function y = squarify(x,d)
% Input x = rectangular matrix
% Output y = square image
% Parameter d indicates interpolation of smaller direction, decimation
% of the larger direction, or compromise of the two
% 's' decimates larger
% 'l' interpolates smaller
% 'm' takes average of both
x = real(squeeze(x));
s = size(x);
if length(s) <3
    if d == 's'
        m = min(s);
        [a, b] = meshgrid(1/s(1):1/s(1):1,1/s(2):1/s(2):1);
        [A, B] = meshgrid(1/m:1/m:1,1/m:1/m:1);
        y = interp2(a,b,x',A,B);
    elseif d == 'l' || ~exist('d')
        m = max(s);
        [a, b] = meshgrid(1/s(1):1/s(1):1,1/s(2):1/s(2):1);
        [A, B] = meshgrid(1/m:1/m:1,1/m:1/m:1);
        y = interp2(a,b,x',A,B);
    elseif d == 'm'
        m = round(mean(s));
        [a, b] = meshgrid(1/s(1):1/s(1):1,1/s(2):1/s(2):1);
        [A, B] = meshgrid(1/m:1/m:1,1/m:1/m:1);
        y = interp2(a,b,x',A,B);
    end
    p = isnan(y(1,:));
    q = isnan(y(:,1));
    z = min([sum(p),sum(q)]);
    y = y(z+1:end,z+1:end); %chops out NaN columns and returns square
else
    if d == 's'
        m = min(s);
        [a, b, c] = meshgrid(1/s(2):1/s(2):1,1/s(1):1/s(1):1,1/s(3):1/s(3):1);
        [A, B, C] = meshgrid(1/m:1/m:1,1/m:1/m:1,1/m:1/m:1);
        y = interp3(a,b,c,x,A,B,C);
    elseif d == 'l' || ~exist('d')
        m = max(s);
        [a, b, c] = meshgrid(1/s(2):1/s(2):1,1/s(1):1/s(1):1,1/s(3):1/s(3):1);
        [A, B, C] = meshgrid(1/m:1/m:1,1/m:1/m:1,1/m:1/m:1);
        y = interp3(a,b,c,x,A,B,C);
    elseif d == 'm'
        m = round(mean(s));
        [a, b, c] = meshgrid(1/s(2):1/s(2):1,1/s(1):1/s(1):1,1/s(3):1/s(3):1);
        [A, B, C] = meshgrid(1/m:1/m:1,1/m:1/m:1,1/m:1/m:1);
        y = interp3(a,b,c,x,A,B,C);
    end
    r = floor(length(y)/2);
    y1 = squeeze(y(:,:,r));
    y2 = squeeze(y(:,r,:));
    y3 = squeeze(y(r,:,:));
    p1 = isnan(y1(1,:));
    q1 = isnan(y1(:,1));
       p2 = isnan(y2(1,:));
    q2 = isnan(y2(:,1));
       p3 = isnan(y3(1,:));
    q3 = isnan(y3(:,1));
    z = [sum(p1) sum(p2) sum(p3) sum(q1) sum(q2) sum(q3)];
    z = z(z<length(y));
    z = max(z);
    y = y(z+1:end,z+1:end,z+1:end); %chops out NaN columns and returns square
end

