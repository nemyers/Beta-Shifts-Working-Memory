function [J,h_theta] = jointcostaccfun(param,x,y,task)

if nargin<4
    task = ones(size(x));
end
sigma  = param(1);
lapses = param(2:3);
if length(param)==4
    mu = param(4);
else
    mu = 0;
end

m = length(y);
h_theta = nan(m,1);

utask = unique(task);
ntask = length(utask);
for itask = 1:ntask
    itrl = task == utask(itask);
    h_theta(itrl) = respfun(x(itrl),mu,sigma,lapses(itask));
end
J = 1/m*sum(-y.*log(h_theta) - (1-y).*log(1-h_theta));

function p = respfun(x,mu,sigma,lapse)
    p = (1-lapse).*normcdf(x,mu,sigma) + lapse/2;
end
end