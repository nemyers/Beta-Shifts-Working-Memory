function  C2 = iem_mahal(data,angs,angspace,foldvar,sigmaind,verbose)
% C2 = iem_mahal(data,angs,angspace,foldvar,sigmaind,verbose)

nTrial   = size(data,1);
nSens    = size(data,2);
nTime    = size(data,3);

if nargin < 6
    verbose = false;
end
if nargin < 5
    sigmaind = ones(1,nTime) > 0;
end
if nargin < 4
    foldvar = 1:nTrial;
end
if nargin < 3
    angspace  = -pi:pi/8:pi-(pi/8);
end

nChan    = length(angspace);
nSigma   = nnz(sigmaind);
C2       = zeros(nChan,nTrial,nTime,'single');
D        = zeros(nChan,nSens,nTime,'single');
folds    = unique(foldvar);
nfold    = length(folds);
X        = [sin(angs) cos(angs) ones(nTrial,1)];

if verbose, ft_progress('init','text'); end
for ifold = 1:nfold    
    if verbose, ft_progress(ifold/nfold); end    
    traintrl  = foldvar ~= folds(ifold);
    testtrl   = find(foldvar == folds(ifold));
    trn_dat   = data(traintrl,:,:);
    Y         = reshape(trn_dat(:,:,sigmaind),[size(trn_dat,1),nSens*nSigma]);    
    res       = Y-X(traintrl,:)*(pinv(X(traintrl,:))*Y);
    res       = reshape(res,[size(trn_dat,1),nSens,nSigma]);
    sigma     = covdiag(mean(res,3));
    
    for trl = testtrl'
        for c = 1:nChan
            ind2 = abs(circ_dist(angs,angs(trl)-angspace(c)))<(pi/8); % bin width of pi/8
            D(c,:,:) = mean(data(traintrl & ind2,:,:),1);
        end

        for time = 1:nTime
            C2(:,trl,time) = pdist2(data(trl,:,time),D(:,:,time),'mahalanobis',sigma);
        end
    end    
end
if verbose, ft_progress('close'); end;
end
