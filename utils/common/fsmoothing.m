function y = fsmoothing(x,b,dsmoothing)

if nargin < 3
    dsmoothing = 2;
end
y = b(1)+diff(b)./(1+exp(-log(1/0.01-1)*2/dsmoothing*x));
