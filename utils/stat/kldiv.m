function KL = kldiv(A,B,normalize,dim)
%KL = kldiv(A,B,normalize,dim)
if nargin < 4
    dim = 2;
end
if nargin < 3
    normalize = false;
end
%%
if normalize
    A = bsxfun(@rdivide,A,sum(A,dim));
    B = bsxfun(@rdivide,B,sum(B,dim));
end

KL = sum(A .* (log(A./B)),dim);
end