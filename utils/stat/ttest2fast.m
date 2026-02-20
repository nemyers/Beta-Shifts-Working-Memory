function [p,tval] = ttest2fast(x,y)
    dim = 1;
    nx = size(x,1);
    ny = size(y,1);
    s2x = nanvar(x,[],dim);
    s2y = nanvar(y,[],dim);
    xmean = nanmean(x,dim);
    ymean = nanmean(y,dim);
    difference = xmean - ymean;
    
    s2xbar = s2x ./ nx;
    s2ybar = s2y ./ ny;
    dfe = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
    se = sqrt(s2xbar + s2ybar);
    se(fix) = 0;
    
    ratio = difference ./ se;
    tval = ratio;
    p = 2 * spm_Tcdf(-abs(ratio),dfe);
end