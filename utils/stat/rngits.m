function [y] = rngits(dim,x,p)
%  RNGITS  Random number generator using inverse transform sampling
%  [y] = rngits(dim,x,p)

if nargin < 3
    error('Missing input argument(s).');
end
if size(x) ~= size(p)
    error('x and p should have the same size.');
end

p = p/trapz(x,p);
c = cumtrapz(x,p);

imin = find(c == 0,1,'last');
if isempty(imin)
    imin = 1;
end
imax = find(c == 1,1,'first');
if isempty(imax)
    imax = length(c);
end
c = c(imin:imax);
x = x(imin:imax);

[c,i] = unique(c);
x = x(i);

y = interp1(c,x,rand(dim));

end