function  [C2,cfg] = mahal_iem_dual_trialwise(dat,theta,cfg)
% [C2,cfg] = mahal_iem_dual_trialwise(dat,theta,cfg)

% authors: M Stokes, E Spaak, N Myers

if nargin < 2
    help mahal_iem_dual_trialwise
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

%flag setting whether covariance is estimated for each time point
if ~isfield(cfg,'covbytime') || isempty(cfg.covbytime)
    cfg.covbytime = false; end
covbytime    = cfg.covbytime;

%boolean vector of time indices to include in the noise covariance estimation
if ~isfield(cfg,'covtime') || isempty(cfg.covtime)
    cfg.covtime = true(1,ntime); end
covtime    = cfg.covtime;

%exclusion mask (opposite of inclusion mask: exclude any matches from
%training set)
if ~isfield(cfg,'exclusion') || isempty(cfg.exclusion)
    cfg.exclusion = [1:ntrl]'; end
exclusion  = cfg.exclusion;

%inclusion mask (only include in training set the trials that match the
%value of the inclusion mask in the test set (assumes each test set has one
%possible inclusion mask value)
if ~isfield(cfg,'inclusion') || isempty(cfg.inclusion)
    cfg.inclusion = true(ntrl,1); end    
inclusion  = cfg.inclusion;

%globabl exclusion mask - set trials to 1 if you never want to include them in training
if ~isfield(cfg,'globalexclusion') || isempty(cfg.globalexclusion)
    cfg.globalexclusion = false(ntrl,1); end
globalexclusion  = cfg.globalexclusion;

%angles (radians)
cfg.theta = theta;
ntheta    = size(theta,2);

%basis set
%create theta basis as cell array so multiple basis sets can be used
if ~isfield(cfg,'thetabasis') || isempty(cfg.thetabasis)
    cfg.thetabasis = {-pi:pi/8:(pi-pi/8)}; end
%backward compatibility with vector inputs
if ~iscell(cfg.thetabasis)
    cfg.thetabasis = {cfg.thetabasis}; end
%revert to vector when using a single basis set for simplicity
if ntheta == 1
    thetabasis = cfg.thetabasis{1};
else
    thetabasis = cfg.thetabasis;
end

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

%recenter tuning curves?
if ~isfield(cfg,'recenter_curve') || isempty(cfg.recenter_curve)
    cfg.recenter_curve = false; end
recenter_curve = cfg.recenter_curve;

%distance function to use
if ~isfield(cfg,'distancefun') || isempty(cfg.distancefun)
    cfg.distancefun = 'mahal'; end
distancefun = cfg.distancefun;

% prepare basis set for design matrix
nChans = [];
for itheta = 1:ntheta
    thetabasis{itheta} = thetabasis{itheta}(:); %make sure the basis set locations are in a column vector
    nChans(1,itheta)   = length(thetabasis{itheta}); %number of virtual channels of output
end
nReg       = sum(nChans); %number of regressors in model (minus constant if used)

%get basis function exponent
if ~isfield(cfg,'basisexponents') || isempty(cfg.basisexponents)
    cfg.basisexponents = nChans-1; end
pows = cfg.basisexponents;

%do test angles match training angles?
if ~isfield(cfg,'testangles') || isempty(cfg.testangles)
    cfg.testangles = theta; end
testangles = cfg.testangles;

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

% set up reverse encoding model parameters (i.e., design matrix)
if add_constant
    constant = 1;
else
    constant = 0;
end

cosfun    = @(theta,mu,pow)((0.5 + 0.5.*cos((theta-mu))).^pow);
C         = zeros(nReg+constant,ntrl);
ireg      = 0;
for itheta = 1:ntheta
    for iChan = 1:nChans(itheta)
        ireg = ireg+1;
        C(ireg,:) = cosfun(theta(:,itheta),thetabasis{itheta}(iChan),pows(itheta));
    end
end
if ~smooth_design_matrix
    iadd = 0;
    for itheta = 1:ntheta
        C([1:nChans(itheta)]+iadd,:) = bsxfun(@minus,C([1:nChans(itheta)]+iadd,:),max(C([1:nChans(itheta)]+iadd,:)))==0;
        iadd = iadd+nChans(itheta);
        %C = (C>0.70) +0;
    end
end
if add_constant
    C(nReg+1,:) = 1;%add constant regressor to model
end
cfg.C = C;
%%
C2      = nan(nReg,ntrl,ntime,nsearch,'single');
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
                    all(bsxfun(@eq,inclusion,inclusion(testtrl(1),:)),2) & ... % keep only trials that match the test trials in the inclusion mask
                    all(bsxfun(@ne,exclusion,exclusion(testtrl(1),:)),2) & ... % keep only trials that differ from test trials in the exclusion mask
                    globalexclusion == 0);                   % exclude some trials from training altogether

    %stratify ?
    
    %get condition means for training data
    pC       = pinv(C(:,traintrl)');
    dtmp     = dat(traintrl,:,:);
    if demean, Dmu = mean(dtmp,1); dtmp = bsxfun(@minus,dtmp,Dmu); end
    dtmp     = reshape(dtmp,[nnz(traintrl) nsens*ntime]);
    b        = pC*dtmp;
    D        = reshape(b,[nReg+constant nsens ntime]);
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
            %if ~covbytime, sigma = covdiag(mean(res(:,:,covtime),3)); end %covdiag or cov
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
    end
end

if recenter_curve
    C    = zeros(nReg,ntrl,ntime,nsearch,'single');
    iadd = 0;
    for itheta = 1:ntheta
        chanlist = 1:nChans(itheta);
        nChan = nChans(itheta);
        for itrl = 1:ntrl
            iang     = testangles(itrl,itheta);
            idst     = circ_dist(thetabasis{itheta},iang);
            [~,imin] = min(abs(idst));
            shiftval = -(ceil(nChan/2) - imin + 1);
            if shiftval<0, shiftval=shiftval+nChan; end
            shiftind     = [chanlist(shiftval+1:nChan) chanlist(1:shiftval)];
            %shiftbasis(itrl,:) = circ_dist(thetabasis(shiftind),iang);
            C(chanlist+iadd,itrl,:) = squeeze(C2(shiftind+iadd,itrl,:,:));
        end
        iadd = iadd+nChans(itheta);
    end
    C2 = C; clear C;
end

if verbose
   %m = sprintf(' - done!\n');
   %fprintf([r m]);
   fprintf([r]);
end