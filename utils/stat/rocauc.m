function [auc,p] = rocauc(x,y,minmax)

if nargin < 3
    minmax = false;
end
if nargin < 2
    error('Missing input argument(s).');
end

auc = [];
p = [];

x = x(:);
y = y(:);

xall = unique(x);
yall = unique(y);
if length(xall) < 2
    auc = 0.5; p = Inf;
    return
end
if length(yall) ~= 2
    error('Invalid input argument(s).');
end

x1 = x(y == yall(2)); n1 = length(x1);
x2 = x(y == yall(1)); n2 = length(x2);

if n1 == 0 || n2 == 0
    return
end

r = tiedrank([x1;x2]);

auc = (sum(r(1:n1))-n1*(n1+1)/2)/(n1*n2);
if minmax
    auc = min(max(auc,1/(n1*n2)),1-1/(n1*n2));
end

if nargout == 2
    p = ranksum(x1,x2);
end

end