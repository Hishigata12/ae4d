function h = blue2(m)
%HOT    Cyan-Blue-Black-Red-Yellow-White color map.
%   HOTCOLD(M) returns an M-by-3 matrix containing a "hotcold" colormap.
%   HOTCOLD, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(hotcold)
%
%   See also HOT, HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, 
%   COLORMAP, RGBPLOT.

%   R. Witte, 8-03-01.
%   Copyright (c) 2001 by Rusty
%   $Revision: 1.0 $  $Date: 2001/11/04 14:33:49 $

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

n = fix(3/8*m);

if mod(m,2)==1
    m = m;
else 
    m = m-1;
end
% b = [(1:n)'/n; ones(m-n,1)];
% g = [zeros(n,1); (1:n)'/n; ones(m-2*n,1).*0.8];
% r = [zeros(2*n,1); (1:m-2*n)'/(m-2*n).*0.8];

r = [zeros(m/2,1); (1:m/2)'/(m/2)*0.8];
b = [(1:m/2)'/(m/2); ones(m/2,1)];
g = [zeros(m/2,1); (1:m/2)'/(m/2)*0.8];

h = [r g b];