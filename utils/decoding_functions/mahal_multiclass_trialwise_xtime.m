function  [C2,cfg] = mahal_twoclass_trialwise_xtime(dat,labels,cfg)
% [C2,cfg] = mahal_twoclass_trialwise_xtime(dat,labels,cfg)
% C2: output, nClass x nTrials x NTrainTime x NTestTime

% authors: M Stokes, E Spaak, N Myers

if nargin < 2
    help mahal_twoclass_trialwise
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

%add constant regressor to design matrix?
if ~isfield(cfg,'add_constant') || isempty(cfg.add_constant)
    cfg.add_constant = false; end
add_constant = cfg.add_constant;

X = dummyvar(labels);
X(:,all(X==repmat(mean(X,1),[size(X,1) 1]),1)) = [];
X = [X ones(length(labels),add_constant)];
nreg = size(X,2);
%%
folds   = unique(foldvar);
folds   = folds(:)';
nfold   = length(folds);

if verbose
   fprintf(' Decoding:'); r = ''; 
end

switch lower(distancefun)
    case {'mahal' 'mahalanobis' 'm' 'mah'}
        for ifold = 1:nfold
            if verbose
               m = sprintf(' Fold %d/%d',ifold,nfold); fprintf([r m]); r = repmat('\b',1,length(m)); 
            end
            testtrl  = find(foldvar == folds(ifold));                % get test trials
            
            traintrl = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                            inclusion == inclusion(testtrl(1)) & ... % keep only trials that match the test trials in the inclusion mask
                            exclusion ~= exclusion(testtrl(1)) & ... % keep only trials that differ from test trials in the exclusion mask
                            globalexclusion == 0);                   % exclude some trials from training altogether

            %get condition means for training data
            pC       = pinv(X(traintrl,:));
            dtmp     = reshape(dat(traintrl,:,:),[nnz(traintrl) nsens*ntime]);
            b        = pC*dtmp;
            res      = dat(traintrl,:,:) - reshape([b'*X(traintrl,:)']',[nnz(traintrl) nsens ntime]);
            D        = reshape(b,[nreg nsens ntime]);

            sigma    = covdiag(mean(res(:,:,covtime),3)); %covdiag or cov

            if add_constant
                Dmu        = D(end,:,:);
                D(end,:,:) = [];
                for ti = 1:ntime
                    for tj = 1:ntime
                        C2(:,testtrl,ti,tj) = pdist2(bsxfun(@minus,dat(testtrl,:,tj),Dmu(1,:,ti)),D(:,:,ti),'mahalanobis',sigma)';
                    end
                end
            else                
                for ti = 1:ntime
                    for tj = 1:ntime
                        C2(:,testtrl,ti,tj) = pdist2(dat(testtrl,:,tj),D(:,:,ti),'mahalanobis',sigma)';
                    end
                end                
            end            
        end
    case {'euclid' 'euclidean' 'e' 'euc'}
        for ifold = 1:nfold
            if verbose
               m = sprintf(' Fold %d/%d',ifold,nfold); fprintf([r m]); r = repmat('\b',1,length(m)); 
            end            
            testtrl  = find(foldvar == folds(ifold));                % get test trials
            
            traintrl = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                            inclusion == inclusion(testtrl(1)) & ... % keep only trials that match the test trials in the inclusion mask
                            exclusion ~= exclusion(testtrl(1)) & ... % keep only trials that differ from test trials in the exclusion mask
                            globalexclusion == 0);                   % exclude some trials from training altogether
            
            %get condition means for training data
            pC       = pinv(X(traintrl,:));
            dtmp     = reshape(dat(traintrl,:,:),[nnz(traintrl) nsens*ntime]);
            D        = reshape(pC*dtmp,[nChan+constant nsens ntime]);

            if add_constant
                Dmu        = D(end,:,:);
                D(end,:,:) = [];
                for ti = 1:ntime
                    for tj = 1:ntime
                        C2(:,testtrl,ti,tj) = pdist2(bsxfun(@minus,dat(testtrl,:,tj),Dmu(1,:,ti)),D(:,:,ti),'euclidean')';
                    end
                end
            else                
                for ti = 1:ntime
                    for tj = 1:ntime
                        C2(:,testtrl,ti,tj) = pdist2(dat(testtrl,:,tj),D(:,:,ti),'euclidean')';
                    end
                end                
            end
        end
end

if verbose
   fprintf(' - done!\n'); 
end