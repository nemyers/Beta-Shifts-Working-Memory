function LL = fit_vonmises(y,x,mu,param)
%LL = fit_vonmises(y,x,mu,param)
%y: vector of responses (0: counterclockwise, 1: clockwise)
%x: angle offset of the probe (in radians!)
%mu: 
sigma = param(1)*ones(size(x));
if length(param) == 2
    lapse = param(2);
    fit_lapse = true;
else
    fit_lapse = false;
end
if length(mu) == 1
    mu = mu*ones(size(x));
end
m     = length(y);

if fit_lapse
    h_theta = respfun_wlapse(x,mu,sigma,lapse);
else
    h_theta = respfun(x,mu,sigma);
end
    

LL = 1/m*sum(-y.*log(h_theta) - (1-y).*log(1-h_theta));

function p = respfun_wlapse(x,mu,sigma,lapse)
         p = (1-lapse).*circ_vmpdf(x,mu,sigma) + lapse/2/pi;  
         %p = (1-lapse).*circ_vmpdf(x,mu,sigma);
end

function p = respfun(x,mu,sigma)
         p = circ_vmpdf(x,mu,sigma);
end

end