ecfunction p = mixpdf(x, mu, K, Pm)
%weighted mixture of von Mises pdf and uniform pdf 

%p = vonmisespdf(x,mu,K)*Pm + 1/(2*pi)*(1-Pm);
p = circ_vmpdf(x,mu,K)*Pm + 1/(2*pi)*(1-Pm);