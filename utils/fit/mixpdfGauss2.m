function p = mixpdfGauss2(x, mu, sigma, Pm, sigma2, Pm2)
%weighted mixture of von Mises pdf and uniform pdf 

p = normpdf(x,mu,sigma)*Pm + normpdf(x,mu,sigma2)*Pm2 + 1/(2*pi)*(1-Pm-Pm2);