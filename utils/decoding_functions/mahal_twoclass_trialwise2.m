function  [C2,cfg] = mahal_twoclass_trialwise2(dat,labels,cfg)
% [C2,cfg] = mahal_twoclass_trialwise2(dat,labels,cfg)

if nargin < 3
    help mahal_twoclass_trialwise2
    error('missing input arguments!')    
end

ntrl  = size(dat,1); % number of timepoints       
nsens = size(dat,2); % number of M/EEG sensors
ntime = size(dat,3);

%set options
if (nargin < 3) || (isempty(cfg))
    cfg = struct;    
end

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

%labels
cfg.labels = labels;

%verbosity
if ~isfield(cfg,'verbose') || isempty(cfg.verbose)
    cfg.verbose = true; end
verbose    = cfg.verbose;

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
%%
C2      = nan(ntrl,ntime,2,'single');
ulabels = unique(labels);
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
    % get same label trials
    traintrl1 = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                    labels == ulabels(1) & ...
                    all(bsxfun(@eq,inclusion,inclusion(testtrl(1),:)),2) & ... % keep only trials that match the test trials in the inclusion mask
                    all(bsxfun(@ne,exclusion,exclusion(testtrl(1),:)),2) & ... % keep only trials that differ from test trials in the exclusion mask
                    globalexclusion == 0);                   % exclude some trials from training altogether

    % get other label trials
    traintrl2 = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                    labels == ulabels(2) & ...
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
    
    %get condition averages
    dtmp     = dat([traintrl1 traintrl2],:,:);
    if demean, Dmu = mean(dtmp,1); dtmp = bsxfun(@minus,dtmp,Dmu); end
    D        = cat(1,mean(dtmp(1:length(traintrl1),:,:),1),mean(dtmp((length(traintrl1)+1):end,:,:),1));
    if add_constant, Dmu = D(end,:,:); D(end,:,:) = []; end
    
    % get test data
    dattest = dat(testtrl,:,:);
    if add_constant || demean
        dattest = bsxfun(@minus,dattest,Dmu); end

    % compute distances to test data
    switch lower(distancefun)
        case {'mahal' 'mahalanobis' 'm' 'mah'}
            % get covariance estimate
            res = bsxfun(@minus,dtmp(1:length(traintrl1),:,:),D(1,:,:));
            res = cat(1,res,bsxfun(@minus,dtmp((length(traintrl1)+1):end,:,:),D(2,:,:)));    
            %if ~covbytime, sigma = covdiag(mean(res(:,:,covtime),3)); end %covdiag or cov                
            if ~covbytime,
                nfeat = size(res,1)*nnz(covtime);
                sigma = covdiag(reshape(permute(res(:,:,covtime),[1 3 2]),nfeat,nsens));
            end
            if covbytime
                for ti = 1:ntime
                    C2(testtrl,ti,:) = pdist2(dattest(:,:,ti),D(:,:,ti),'mahalanobis',covdiag(res(:,:,ti)));
                end                    
            else
                for ti = 1:ntime
                    C2(testtrl,ti,:) = pdist2(dattest(:,:,ti),D(:,:,ti),'mahalanobis',sigma);
                end
            end

        case {'euclid' 'euclidean' 'e' 'euc'}
            for ti = 1:ntime                  
                C2(testtrl,ti,:) = pdist2(dattest(:,:,ti),D(:,:,ti),'euclidean');
            end
    end
end
cfg.muvec = muvec;

C2(labels == ulabels(2),:,:) = C2(labels == ulabels(2),:,[2 1]);

if verbose
   fprintf(' - done!\n'); 
end