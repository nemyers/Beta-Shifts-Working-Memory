function [betas,residuals] = massGLM(X,y,add_offset,rank_regression)
%function [betas,residuals] = massGLM(X,y,add_offset,rank_regression)
%efficiently calculates betas for mass univariate GLM
%
%inputs: X (design matrix) should be nTrials*nRegressors
%        y (data) should be nTrials*nChannels*nTimepoints(*nFrequencies)
%        add_offset should be true to add a constant regressor to
%        the design matrix, false if not (default = true)
%        rank_regression should be true to rank order regressors and data
%        before regression (default = false)
%
%outputs: betas is nRegressors*nChannels*nTimepoints(*nFrequencies)

if nargin < 4
    rank_regression = false;
end
if nargin < 3
    add_offset = true;
end

%make sure inputs are in order
if size(X,2)>size(X,1)
    X = X';
    warning('massGLM: transposing design matrix')
end
if any(all(X==1,1)) & add_offset
    add_offset = false;
    warning('massGLM: design matrix already appears to have a constant regressor.')
end

%get data and design matrix dimensions
dims  = size(y);
ntrl  = size(X,1);
if dims(1) ~= ntrl
    error('massGLM: trial numbers in X and y don''t match!');
end
X     = [ones(ntrl,add_offset) X];
nregs = size(X,2);

%reshape data matrix so it is nTrials*nFeatures
%(=channelsXtimepointsXfrequencies...)
y  = reshape(y,[dims(1) prod(dims(2:end))]);

%rank order X and y?
if rank_regression
    [~,I] = sort(y);
    [~,y] = sort(I);
    [~,I] = sort(X);
    [~,X] = sort(I);
end

%estimate betas
pX    = pinv(X);
betas = pX*y;

%get residuals if requested
if nargout==2
    residuals = y-X*betas;
    residuals = reshape(residuals,dims);
end

%reshape betas to sensible output
betas = reshape(betas,[nregs dims(2:end)]);

end