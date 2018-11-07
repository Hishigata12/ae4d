function h = hotcold(m)
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

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n = fix(1/4*m);p=fix(1/4*m);

r = [zeros(2*n,1);0;(1:n)'/n; ones(n,1)];
g = [1-(0:n-1)'/n;zeros(2*n,1);0;(1:n)'/n];
b = [ones(n,1);1-(0:n-1)'/n;0;zeros(2*n,1)];

h = [r g b];
