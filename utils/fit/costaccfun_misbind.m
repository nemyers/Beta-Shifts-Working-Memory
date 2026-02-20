function J = costaccfun_misbind(param,x,y,usepriors)

if nargin < 4
    usepriors = false;
end

mu    = 0;
sigma = param(1);
lapse = param(2);
msbnd = param(3);
if length(param)==4
    sigma2 = param(4);
else
    sigma2 = sigma;
end

if usepriors
    %for Lapse rate: Beta distribution
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
    p = (1-lapse/2-msbnd).*normcdf(x(:,1),mu,sigma) + msbnd.*normcdf(x(:,2),mu,sigma2) + lapse/2;
end
end