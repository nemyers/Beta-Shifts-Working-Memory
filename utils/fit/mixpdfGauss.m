function p = mixpdfGauss2(x, mu, sigma, Pm)
%weighted mixture of von Mises pdf and uniform pdf 

p = normpdf(x,mu,sigma)*Pm + 1/(2*pi)*(1-Pm);