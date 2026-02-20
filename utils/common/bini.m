function [ibin] = bini(x,nbin,pbin)

if nargin < 3
    pbin = 1/nbin;
end
if nargin < 2
    error('Missing input arguments!');
end

if ~isvector(x)
    error('Wrong input size!');
end

xbinbeg = quantile(x(:),linspace(0,1-pbin,nbin));
xbinend = quantile(x(:),linspace(pbin,1,nbin));

ibin = false(numel(x),nbin);
for i = 1:nbin
    ibin(:,i) = x(:) >= xbinbeg(i) & x(:) <= xbinend(i);
end

end