function  [out_proj,variance_proj,cfg] = mante_proj(dat,cfg)
% [out_proj,variance_proj,cfg] = mante_proj(dat,cfg)

if nargin < 2
    help mante_proj
    error('missing input arguments!')    
end

ntrl  = size(dat,1); % number of timepoints       
nsens = size(dat,2); % number of M/EEG sensors
ntime = size(dat,3);

%set options
if (nargin < 2) || (isempty(cfg))
    cfg = struct;    
end

%folding variable
if ~isfield(cfg,'foldvar') || isempty(cfg.foldvar)
    cfg.foldvar = [1:ntrl]'; end
foldvar    = cfg.foldvar;

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

%add constant regressor to design matrix?
if ~isfield(cfg,'add_constant') || isempty(cfg.add_constant)
    cfg.add_constant = false; end
add_constant = cfg.add_constant;

%run PCA?
if ~isfield(cfg,'run_pca') || isempty(cfg.run_pca)
    cfg.run_pca = true; end
run_pca = cfg.run_pca;

%PCA dimensions
if ~isfield(cfg,'npc') || isempty(cfg.npc)
    cfg.npc = 12; end
npc = cfg.npc;

%orthogonalize betas?
if ~isfield(cfg,'orthogonalize_betas') || isempty(cfg.orthogonalize_betas)
    cfg.orthogonalize_betas = true; end
orthogonalize_betas = cfg.orthogonalize_betas;
    
% set up reverse encoding model parameters (i.e., design matrix)
C = cfg.C;
% prepare basis set for design matrix
nChan      = size(C,1); %number of virtual channels of output

if add_constant
    C(nChan+1,:) = 1;%add constant regressor to model
end
cfg.C = C;

%make design matrix for GLM
X  = cfg.X; %design matrix for circular-linear multiple regression
nX = size(X,2);
%%
if verbose
   fprintf(' Running Mante:'); r = ''; 
end
%pre-allocate
folds   = unique(foldvar);
folds   = folds(:)';
nfold   = length(folds);
%run pca if desired
if run_pca
    pcdat = zeros(ntrl,npc,ntime,'single');
    [coeff,~,~,~,expl_var] = pca(reshape(permute(dat,[2 3 1]),[nsens ntime*ntrl])','Algorithm','svd');
    for it = 1:ntime %project each trial into PC subspace
        pcdat(:,:,it) = dat(:,:,it)*coeff(:,1:npc);
    end
    dat = pcdat; clear pcdat
    nsens = npc;
else
    expl_var = NaN;
end
%loop over x-validation folds
for ifold = 1:nfold
    if verbose
       m = sprintf(' Fold %d/%d',ifold,nfold); fprintf([r m]); r = repmat('\b',1,length(m)); 
    end
    bintesttrl = foldvar == folds(ifold);
    testtrl    = find(bintesttrl)';                           % get test trials
    % get same task trials
    traintrl  = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                    all(bsxfun(@eq,inclusion,inclusion(testtrl(1),:)),2) & ... % keep only trials that match the test trials in the inclusion mask
                    all(bsxfun(@ne,exclusion,exclusion(testtrl(1),:)),2) & ... % keep only trials that differ from test trials in the exclusion mask
                    globalexclusion == 0);                    % exclude some trials from training altogether

    testtrl = testtrl';

    ntrain = length(traintrl);
    ntest  = length(testtrl);

    % get betas
    nvox    = nsens*ntime;
    bdat    = reshape(dat(traintrl,:,:),[ntrain nvox]); %reorder data matrix into 2D ( trials * (sensors*timepoints) ) for GLM
    bdatmu  = mean(bdat,1);
    bdatsd  = std(bsxfun(@minus,bdat,bdatmu),[],1);
    zbdat   = bsxfun(@rdivide,bsxfun(@minus,bdat,bdatmu),bdatsd);

    pX      = pinv(X(traintrl,:));
    nreg    = size(pX,1);
    betacur = single(pX*zbdat); %solve GLM
    betacur = reshape(betacur,[nreg nsens ntime]);

    bdatmu  = reshape(bdatmu,[1 nsens ntime]);
    bdatsd  = reshape(bdatsd,[1 nsens ntime]);

    % orthogonalize betas
    if orthogonalize_betas & nreg>1
        keep_magnitude = true;
        for it = 1:ntime
            betacur(:,:,it) = symmetric_orthogonalise(betacur(:,:,it)',keep_magnitude)';
        end
    end

    % get means from test data
    muerf  = NaN(nsens,nChan,ntime,'single');
    NotNaN = true(nChan,1);
    %mcur  = mean(dat(bintesttrl,:,:),1);
    %sdcur = std(dat(bintesttrl,:,:),[],1);
    for iang = 1:nChan
        iitrl           = C(iang,:)'==1 & bintesttrl;
        NotNaN(iang,1)  = nnz(iitrl)>0;
        if NotNaN(iang)
            %curdat          = bsxfun(@rdivide,bsxfun(@minus,dat(iitrl,:,:),mcur),sdcur);                    
            curdat          = bsxfun(@rdivide,bsxfun(@minus,dat(iitrl,:,:),bdatmu),bdatsd);
            muerf(:,iang,:) = mean(curdat,1);
        end
    end

    % project test data into task subspace from training data
    for it = 1:ntime
        curproj(1:nreg,:,it,ifold) = betacur(:,:,it)*muerf(:,:,it);

        %get eigenvalues of covariance matrix of full data and
        %subspace projection
        stimind = [ones(16,1);zeros(16,1)];
        [~,d]      = eig(cov(muerf(:,stimind & NotNaN,it)'));
        fullvar    = sum(diag(d));

        [~,d]      = eig(cov(curproj(1:2,stimind & NotNaN,it,ifold)'));
        reducedvar = sum(diag(d));

        [~,d]      = eig(cov(curproj(3:4,stimind & NotNaN,it,ifold)'));
        reducedvar2 = sum(diag(d));

        cur_variance_proj(:,it,ifold,1) = [fullvar reducedvar reducedvar2];

        stimind = [zeros(16,1);ones(16,1)];
        [~,d]      = eig(cov(muerf(:,stimind & NotNaN,it)'));
        fullvar    = sum(diag(d));

        [~,d]      = eig(cov(curproj(1:2,stimind & NotNaN,it,ifold)'));
        reducedvar = sum(diag(d));

        [~,d]      = eig(cov(curproj(3:4,stimind & NotNaN,it,ifold)'));
        reducedvar2 = sum(diag(d));

        cur_variance_proj(:,it,ifold,2) = [fullvar reducedvar reducedvar2];
    end
end

cfg.expl_var  = expl_var;
out_proj    = curproj;
variance_proj = cur_variance_proj;


if verbose
   fprintf(' - done!\n'); 
end