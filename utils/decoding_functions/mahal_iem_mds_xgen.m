function  [distmat,cfg] = mahal_iem_mds_xgen(dat,theta,task,cfg)
% [distmat,cfg] = mahal_iem_mds_xgen(dat,theta,task,cfg)

if nargin < 3
    help mahal_iem_mds_xgen
    error('missing input arguments!')    
end

ntrl  = size(dat,1); % number of timepoints       
nsens = size(dat,2); % number of M/EEG sensors
ntime = size(dat,3);

%set options
if (nargin < 3) || (isempty(cfg))
    cfg = struct;    
end

%basis set
if ~isfield(cfg,'thetabasis') || isempty(cfg.thetabasis)
    cfg.thetabasis = -pi:pi/8:(pi-pi/8); end
thetabasis = cfg.thetabasis;

%folding variable
if ~isfield(cfg,'foldvar') || isempty(cfg.foldvar)
    cfg.foldvar = [1:ntrl]'; end
foldvar    = cfg.foldvar;

%flag setting whether covariance is estimated for each time point
if ~isfield(cfg,'covbytime') || isempty(cfg.covbytime)
    cfg.covbytime = false; end
covbytime    = cfg.covbytime;

%boolean vector of time indices to include in the noise covariance estimation
%(only relevant if covbytime = false)
if ~isfield(cfg,'covtime') || isempty(cfg.covtime)
    cfg.covtime = true(1,ntime); end
covtime    = cfg.covtime;

%inclusion mask (only include in training set the trials that match the
%value of the inclusion mask in the test set (assumes each test set has one
%possible inclusion mask value)
if ~isfield(cfg,'inclusion') || isempty(cfg.inclusion)
    cfg.inclusion = true(ntrl,1); end    
inclusion  = cfg.inclusion;

%exclusion mask (opposite of inclusion mask: exclude any matches from
%training set)
if ~isfield(cfg,'exclusion') || isempty(cfg.exclusion)
    cfg.exclusion = [1:ntrl]'; end
exclusion  = cfg.exclusion;

%globabl exclusion mask - set trials to 1 if you never want to include them in training
if ~isfield(cfg,'globalexclusion') || isempty(cfg.globalexclusion)
    cfg.globalexclusion = false(ntrl,1); end
globalexclusion  = cfg.globalexclusion;

%angles (radians)
cfg.theta = theta;

%verbosity
if ~isfield(cfg,'verbose') || isempty(cfg.verbose)
    cfg.verbose = true; end
verbose    = cfg.verbose;

%smooth design matrix?
if ~isfield(cfg,'smooth_design_matrix') || isempty(cfg.smooth_design_matrix)
    cfg.smooth_design_matrix = true; end
smooth_design_matrix = cfg.smooth_design_matrix;

%demean data? (instead of adding a constant regressor)
if ~isfield(cfg,'demean') || isempty(cfg.demean)
    cfg.demean = false; end
demean = cfg.demean;
if demean, cfg.add_constant = false; end

%add constant regressor to design matrix?
if ~isfield(cfg,'add_constant') || isempty(cfg.add_constant)
    cfg.add_constant = true; end
add_constant = cfg.add_constant;

%distance function to use
if ~isfield(cfg,'distancefun') || isempty(cfg.distancefun)
    cfg.distancefun = 'mahal'; end
distancefun = cfg.distancefun;

% prepare basis set for design matrix
thetabasis = thetabasis(:); %make sure the basis set locations are in a column vector
nChan      = length(thetabasis); %number of virtual channels of output

%get basis function exponent
if ~isfield(cfg,'basisexponent') || isempty(cfg.basisexponent)
    cfg.basisexponent = nChan-1; end
pow = cfg.basisexponent;
    
% set up reverse encoding model parameters (i.e., design matrix)
if add_constant
    constant = 1;
else
    constant = 0;
end
cosfun    = @(theta,mu,pow)((0.5 + 0.5.*cos((theta-mu))).^pow);
C         = zeros(nChan+constant,ntrl);
for iChan = 1:nChan
    C(iChan,:) = cosfun(theta,thetabasis(iChan),pow);
end
if ~smooth_design_matrix
    C(1:nChan,:) = bsxfun(@minus,C(1:nChan,:),max(C(1:nChan,:)))==0;
    %C = (C>0.70) +0;
end
if add_constant
    C(nChan+1,:) = 1;%add constant regressor to model
end
cfg.C = C;
%%
folds   = unique(foldvar);
folds   = folds(:)';
nfold   = length(folds);
muvec   = [];

if verbose
   fprintf(' Decoding:'); r = ''; 
end

for ifold = 1:nfold
    if verbose
       m = sprintf(' Fold %d/%d',ifold,nfold); fprintf([r m]); r = repmat('\b',1,length(m)); 
    end
    testtrl  = find(foldvar == folds(ifold))';                % get test trials
    % get same task trials
    traintrl1 = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                    task == task(testtrl(1)) & ...
                    all(bsxfun(@eq,inclusion,inclusion(testtrl(1),:)),2) & ... % keep only trials that match the test trials in the inclusion mask
                    all(bsxfun(@ne,exclusion,exclusion(testtrl(1),:)),2) & ... % keep only trials that differ from test trials in the exclusion mask
                    globalexclusion == 0);                   % exclude some trials from training altogether

    % get other task trials
    traintrl2 = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                    task ~= task(testtrl(1)) & ...
                    all(bsxfun(@eq,inclusion,inclusion(testtrl(1),:)),2) & ... % keep only trials that match the test trials in the inclusion mask
                    all(bsxfun(@ne,exclusion,exclusion(testtrl(1),:)),2) & ... % keep only trials that differ from test trials in the exclusion mask
                    globalexclusion == 0);                   % exclude some trials from training altogether

    % match number of training trials
    trllen    = [length(traintrl1) length(traintrl2)];
    mintrl    = min(trllen);
    truncind  = find(trllen>mintrl);
    if truncind == 1
        itrl1     = randperm(trllen(1));
        traintrl1 = traintrl1(itrl1(1:mintrl));
        traintrl1 = sort(traintrl1);
    elseif truncind == 2
        itrl2     = randperm(trllen(2));
        traintrl2 = traintrl2(itrl2(1:mintrl));
        traintrl2 = sort(traintrl2);
    end
    
    trldiff1 = mean(abs(bsxfun(@minus,[traintrl1],testtrl)),2);
    mu(1)    = mean(trldiff1,1);
    trldiff2 = mean(abs(bsxfun(@minus,[traintrl2],testtrl)),2);
    mu(2)    = mean(trldiff2,1);
    
    maxiter      = 1e3;
    diffthresh   = 5;
    iter         = 1;
    maxreduction = 0.20; %maximum proportion of trials to throw out
    while abs(diff(mu))>diffthresh
        if iter>maxiter || (length(traintrl1)./mintrl < (1-maxreduction))
            break;
        end
        %which mean is larger?
        musign = sign(diff(mu));
        %remove trial pair
        if musign == -1 %traintrl1 is larger
            [~,rmind] = max(trldiff1);
            traintrl1(rmind) = [];
            [~,rmind] = min(trldiff2);
            traintrl2(rmind) = [];
        elseif musign == +1 %traintrl2 is larger
            [~,rmind] = min(trldiff1);
            traintrl1(rmind) = [];
            [~,rmind] = max(trldiff2);
            traintrl2(rmind) = [];
        elseif musign == 0 %matched!
            break;
        end
        %update mu
        trldiff1 = mean(abs(bsxfun(@minus,[traintrl1],testtrl)),2);
        mu(1)    = mean(trldiff1,1);
        trldiff2 = mean(abs(bsxfun(@minus,[traintrl2],testtrl)),2);
        mu(2)    = mean(trldiff2,1);
        iter = iter+1;
    end
    muvec(ifold,:) = mu;
    testtrl = testtrl';
    
    %loop over tasks
    for itask = 1:2
        if itask == 1, traintrl = traintrl1; else traintrl = traintrl2; end
        
        % get condition means for training data
        pC       = pinv(C(:,traintrl)');
        dtmp     = dat(traintrl,:,:);
        if demean, Dmu = mean(dtmp,1); dtmp = bsxfun(@minus,dtmp,Dmu); end
        dtmp     = reshape(dtmp,[nnz(traintrl) nsens*ntime]);
        b        = pC*dtmp;
        D        = reshape(b,[nChan+constant nsens ntime]);
        if add_constant, Dmu = D(end,:,:); D(end,:,:) = []; end
        
        % get test data
        dattest = dat(testtrl,:,:);
        if add_constant || demean
            dattest = bsxfun(@minus,dattest,Dmu); end
        
        % compute distances to test data
        switch lower(distancefun)
            case {'mahal' 'mahalanobis' 'm' 'mah'}
                % get covariance estimate
                res      = reshape(dtmp - [b'*C(:,traintrl)]',[nnz(traintrl) nsens ntime]);
                if ~covbytime, sigma = covdiag(mean(res(:,:,covtime),3)); end %covdiag or cov                
       
                if covbytime
                    for ti = 1:ntime
                        dists = pdist(D(:,:,ti),'mahalanobis',covdiag(res(:,:,ti)))';
                        distmat(:,:,ti,itask,ifold) = mdscale(squareform(dists),2);                        
                    end                    
                else
                    for ti = 1:ntime
                        dists = pdist(D(:,:,ti),'mahalanobis',sigma)';
                        distmat(:,:,ti,itask,ifold) = mdscale(squareform(dists),2);
                    end
                end

            case {'euclid' 'euclidean' 'e' 'euc'}
                for ti = 1:ntime                  
                    dists = pdist(D(:,:,ti),'euclidean',sigma)';
                    distmat(:,:,ti,itask,ifold) = mdscale(squareform(dists),2);
                end
        end                        
    end
end
cfg.muvec = muvec;
distmat   = mean(distmat,5);

if verbose
   fprintf(' - done!\n'); 
end