function [J,h_theta] = costaccfun(param,x,y,usepriors)

if nargin < 4
    usepriors = false;
end

sigma = param(1);
lapse = param(2);
if length(param)==3
    mu = param(3);
else
    mu = 0;
end

if usepriors
    %for Lapse rate: Beta distributio
    bp         = [1.1 1.1]; %uninformative prior
    priorLapse = betapdf(lapse*1,bp(1),bp(2));
    
    gp         = [1.1 50];
    priorSigma = gampdf(sigma,gp(1),gp(2));
    priorSigma = 1;
else
    priorSigma = 1;
    priorLapse = 1;
end

m     = length(y);

h_theta = respfun(x,mu,sigma,lapse);

J = 1/m*sum(-y.*log(h_theta) - (1-y).*log(1-h_theta)) - log(priorSigma) - log(priorLapse);
%J = 1/m*sum(log(h_theta)) - log(priorSigma) - log(priorLapse);


function p = respfun(x,mu,sigma,lapse)
    p = (1-lapse).*normcdf(x,mu,sigma) + lapse/2;
end
end