function  [C2,cfg] = mahal_multiclass_trialwise(dat,labels,cfg)
% [C2,cfg] = mahal_multiclass_trialwise(dat,labels,cfg)

% authors: M Stokes, E Spaak, N Myers

if nargin < 2
    help mahal_multiclass_trialwise
    error('missing input arguments!')    
end

ntrl  = size(dat,1); % number of timepoints       
nsens = size(dat,2); % number of M/EEG sensors
ntime = size(dat,3);

%set options
if (nargin < 3) || (isempty(cfg))
    cfg = struct;    
end

%fold variable
if ~isfield(cfg,'foldvar') || isempty(cfg.foldvar)
    cfg.foldvar = [1:ntrl]'; end
foldvar    = cfg.foldvar;

%boolean vector of time indices to include in the noise covariance estimation
if ~isfield(cfg,'covtime') || isempty(cfg.covtime)
    cfg.covtime = true(1,ntime); end
covtime    = cfg.covtime;

%flag setting whether covariance is estimated for each time point
if ~isfield(cfg,'covbytime') || isempty(cfg.covbytime)
    cfg.covbytime = false; end
covbytime    = cfg.covbytime;

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

%verbosity
if ~isfield(cfg,'verbose') || isempty(cfg.verbose)
    cfg.verbose = true; end
verbose    = cfg.verbose;

%distance function to use
if ~isfield(cfg,'distancefun') || isempty(cfg.distancefun)
    cfg.distancefun = 'mahal'; end
distancefun = cfg.distancefun;

%run searchlight?
if ~isfield(cfg,'chanlist') || isempty(cfg.chanlist)
    cfg.chanlist = {1:nsens}; end
chanlist = cfg.chanlist;
nsearch  = length(chanlist);
do_searchlight = nsearch>1;

%use shrinkage to estimate covariance?
if ~isfield(cfg,'covshrinkage') || isempty(cfg.covshrinkage)
    cfg.covshrinkage = true; end
%if do_searchlight, cfg.covshrinkage = false; end
covshrinkage = cfg.covshrinkage;

%add constant regressor to design matrix?
if ~isfield(cfg,'add_constant') || isempty(cfg.add_constant)
    cfg.add_constant = false; end
add_constant = cfg.add_constant;

%demean data? (instead of adding a constant regressor)
if ~isfield(cfg,'demean') || isempty(cfg.demean)
    cfg.demean = false; end
demean = cfg.demean;
if demean, cfg.add_constant = false; end

%match number of trials per training label
if ~isfield(cfg,'match_training_label_nums') || isempty(cfg.match_training_label_nums)
    cfg.match_training_label_nums = false; end
match_training_label_nums = cfg.match_training_label_nums;

%set min number of trials per training label
if ~isfield(cfg,'min_training_label_nums') || isempty(cfg.min_training_label_nums)
    cfg.min_training_label_nums = 0;
    if match_training_label_nums==1,
        cfg.min_training_label_nums = 1; end
end
min_training_label_nums = cfg.min_training_label_nums;

X = dummyvar(labels);
X(:,all(X==repmat(mean(X,1),[size(X,1) 1]),1)) = [];
X = [X ones(length(labels),add_constant)];
nreg = size(X,2);
nlab = nreg-add_constant;
ulab = 1:nlab;
%%
folds   = unique(foldvar);
folds   = folds(:)';
nfold   = length(folds);

if verbose
   r = ''; 
end

for ifold = 1:nfold
    if verbose
       m = sprintf(' Decoding: Fold %d/%d',ifold,nfold); fprintf([r m]); r = repmat('\b',1,length(m)); 
    end
    testtrl  = find(foldvar == folds(ifold));                % get test trials

    traintrl = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                    inclusion == inclusion(testtrl(1)) & ... % keep only trials that match the test trials in the inclusion mask
                    exclusion ~= exclusion(testtrl(1)) & ... % keep only trials that differ from test trials in the exclusion mask
                    globalexclusion == 0);                   % exclude some trials from training altogether

    trainlabels = labels(traintrl);
    ntrain = hist(trainlabels,ulab);
    mintrl = min(ntrain);            
    if mintrl>=min_training_label_nums
        if match_training_label_nums
            %get number of training trials per label                
            imin   = find(ntrain==mintrl);
            tmptrain = [];
            for ilab = 1:nlab
                ii = find(trainlabels == ulab(ilab));
                currtrl = traintrl(ii);
                if ilab==imin
                    tmptrain = [tmptrain;currtrl];
                else
                    currtrl = currtrl(randperm(length(currtrl)));
                    tmptrain = [tmptrain;currtrl(1:mintrl)];
                end                
            end
            traintrl = tmptrain;
        end

        %get condition means for training data
        pC       = pinv(X(traintrl,:));
        dtmp     = dat(traintrl,:,:);
        if demean, Dmu = mean(dtmp,1); dtmp = bsxfun(@minus,dtmp,Dmu); end
        dtmp     = reshape(dtmp,[nnz(traintrl) nsens*ntime]);
        b        = pC*dtmp;
        D        = reshape(b,[nreg nsens ntime]);
        if add_constant, Dmu = D(end,:,:); D(end,:,:) = []; end

        % get test data
        dattest = dat(testtrl,:,:);
        if add_constant || demean
            dattest = bsxfun(@minus,dattest,Dmu); end

        switch lower(distancefun)
            case {'mahal' 'mahalanobis' 'm' 'mah'}
                res      = reshape(dtmp - [b'*X(traintrl,:)']',[nnz(traintrl) nsens ntime]);

                if ~covbytime, nfeat = nnz(traintrl)*nnz(covtime);
                    if covshrinkage, sigma = covdiag(reshape(permute(res(:,:,covtime),[1 3 2]),nfeat,nsens)); 
                    else sigma = cov(reshape(permute(res(:,:,covtime),[1 3 2]),nfeat,nsens)); end; end

                if covbytime
                    for ti = 1:ntime
                        if covshrinkage, sigma = covdiag(res(:,:,ti)); else, sigma = cov(res(:,:,ti)); end
                        if do_searchlight
                            for isearch = 1:nsearch
                                ichan = chanlist{isearch};
                                C2(:,testtrl,ti,isearch) = pdist2(dattest(:,ichan,ti),D(:,ichan,ti),'mahalanobis',sigma(ichan,ichan))';
                            end
                        else
                                C2(:,testtrl,ti,1) = pdist2(dattest(:,:,ti),D(:,:,ti),'mahalanobis',sigma)';
                        end
                    end                    
                else
                    for ti = 1:ntime
                        if do_searchlight
                            for isearch = 1:nsearch
                                ichan = chanlist{isearch};
                                sigtmp = sigma(ichan,ichan);
                                C2(:,testtrl,ti,isearch) = pdist2(dattest(:,ichan,ti),D(:,ichan,ti),'mahalanobis',sigtmp)';
                            end
                        else
                                C2(:,testtrl,ti,1) = pdist2(dattest(:,:,ti),D(:,:,ti),'mahalanobis',sigma)';
                        end
                    end
                end
            case {'euclid' 'euclidean' 'e' 'euc'}
                for ti = 1:ntime
                    if do_searchlight
                        for isearch = 1:nsearch
                            ichan = chanlist{isearch};
                            C2(:,testtrl,ti,isearch) = pdist2(dattest(:,ichan,ti),D(:,ichan,ti),'euclidean')';
                        end
                    else
                            C2(:,testtrl,ti,1) = pdist2(dattest(:,:,ti),D(:,:,ti),'euclidean')';
                    end
                end
            case {'lda' 'lineardiscriminant' 'l'}
                res      = reshape(dtmp - [b'*X(traintrl,:)']',[nnz(traintrl) nsens ntime]);

                if ~covbytime, nfeat = nnz(traintrl)*nnz(covtime);
                    if covshrinkage, sigma = covdiag(reshape(permute(res(:,:,covtime),[1 3 2]),nfeat,nsens)); 
                    else sigma = cov(reshape(permute(res(:,:,covtime),[1 3 2]),nfeat,nsens)); end; end

                if covbytime
                    for ti = 1:ntime
                        if covshrinkage, sigma = covdiag(res(:,:,ti)); else, sigma = cov(res(:,:,ti)); end
                        if do_searchlight
                            for isearch = 1:nsearch
                                ichan = chanlist{isearch};
                                C2(:,testtrl,ti,isearch) = pdist2(dattest(:,ichan,ti),D(:,ichan,ti),'mahalanobis',sigma(ichan,ichan))';
                            end
                        else
                                C2(:,testtrl,ti,1) = pdist2(dattest(:,:,ti),D(:,:,ti),'mahalanobis',sigma)';
                        end
                    end                    
                else
                    for ti = 1:ntime
                        if do_searchlight
                            for isearch = 1:nsearch
                                ichan = chanlist{isearch};
                                sigtmp = sigma(ichan,ichan);
                                C2(:,testtrl,ti,isearch) = pdist2(dattest(:,ichan,ti),D(:,ichan,ti),'mahalanobis',sigtmp)';
                            end
                        else
                                C2(:,testtrl,ti,1) = pdist2(dattest(:,:,ti),D(:,:,ti),'mahalanobis',sigma)';
                        end
                    end
                end
        end
    else
        C2(:,testtrl,1:ntime,1:nsearch) = NaN;
    end
end
if verbose
   fprintf([r ]);
end