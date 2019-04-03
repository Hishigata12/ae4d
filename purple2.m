function h = purple2(m)
if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

n = fix(3/8*m);

if mod(m,2) == 0
    m = m;
else 
    m = m-1;
end


% r = [(1:n)'/n; ones(m-n,1)];
% g = [zeros(n,1); (1:n)'/n; ones(m-2*n,1)];
% b = [zeros(2*n,1); (1:m-2*n)'/(m-2*n)];
r = (1:m)'/m*0.8;
b = (1:m)'/m*0.8;
g = [zeros(m/2,1); (1:m/2)'/(m/2)*0.8];
h = [r g b];