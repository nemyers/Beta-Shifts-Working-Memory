function  [variance_proj,cfg] = pca_proj_xgen(dat,theta,task,cfg)
% [variance_proj,cfg] = pca_proj_xgen(dat,theta,task,cfg)

if nargin < 3
    help pca_proj_xgen
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
    cfg.smooth_design_matrix = false; end
smooth_design_matrix = cfg.smooth_design_matrix;

%add constant regressor to design matrix?
if ~isfield(cfg,'add_constant') || isempty(cfg.add_constant)
    cfg.add_constant = false; end
add_constant = cfg.add_constant;

% prepare basis set for design matrix
thetabasis = thetabasis(:); %make sure the basis set locations are in a column vector
nChan      = length(thetabasis); %number of virtual channels of output

%get basis function exponent
if ~isfield(cfg,'basisexponent') || isempty(cfg.basisexponent)
    cfg.basisexponent = nChan-1; end
pow = cfg.basisexponent;

%return only design matrix and quit?
if ~isfield(cfg,'returnConly') || isempty(cfg.returnConly)
    cfg.returnConly = false; end
returnConly = cfg.returnConly;
    
%match distributions of training trials in terms of distance to test
%trials?
if ~isfield(cfg,'match_train_trial_distance_to_test_trial') || isempty(cfg.match_train_trial_distance_to_test_trial)
    cfg.match_train_trial_distance_to_test_trial = true; end
match_train_trial_distance_to_test_trial = cfg.match_train_trial_distance_to_test_trial;

%match numbers of training trials?
if ~isfield(cfg,'match_training_trial_numbers') || isempty(cfg.match_training_trial_numbers)
    cfg.match_training_trial_numbers = true; end
match_training_trial_numbers = cfg.match_training_trial_numbers;

%run PCA?
if ~isfield(cfg,'run_pca') || isempty(cfg.run_pca)
    cfg.run_pca = true; end
run_pca = cfg.run_pca;

%PCA dimensions
if ~isfield(cfg,'npc') || isempty(cfg.npc)
    cfg.npc = 12; end
npc = cfg.npc;
    
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

%make design matrix for GLM
if ~isfield(cfg,'X') || isempty(cfg.X)
    X  = [sin(theta) cos(theta) task]; %design matrix for circular-linear multiple regression
    cfg.X = X; end
X  = cfg.X;
nX = size(X,2);
%%
if returnConly
    variance_proj = [];
    cfg.expl_var = [];
else
    if verbose
       fprintf(' Running Mante:'); r = ''; 
    end
    %pre-allocate
    folds   = unique(foldvar);
    folds   = folds(:)';
    nfold   = length(folds);
    muvec   = [];
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
        testtask   = task(testtrl(1));
        % get same task trials
        traintrl1 = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                        task == testtask & ...
                        all(bsxfun(@eq,inclusion,inclusion(testtrl(1),:)),2) & ... % keep only trials that match the test trials in the inclusion mask
                        all(bsxfun(@ne,exclusion,exclusion(testtrl(1),:)),2) & ... % keep only trials that differ from test trials in the exclusion mask
                        globalexclusion == 0);                    % exclude some trials from training altogether

        % get other task trials
        traintrl2 = find(foldvar ~= folds(ifold) & ...            % use other trials for training
                        task ~= testtask & ...
                        all(bsxfun(@eq,inclusion,inclusion(testtrl(1),:)),2) & ... % keep only trials that match the test trials in the inclusion mask
                        all(bsxfun(@ne,exclusion,exclusion(testtrl(1),:)),2) & ... % keep only trials that differ from test trials in the exclusion mask
                        globalexclusion == 0);                    % exclude some trials from training altogether

        % match number of training trials
        if match_training_trial_numbers
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
        end

        % match distance of training trials to test trials
        if match_train_trial_distance_to_test_trial
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
        end
        trainall = [traintrl1;traintrl2];
        
        testtrl = testtrl';

        %loop over tasks
        for itask = 1:3
            if     itask == 1, traintrl = traintrl1; %same task 
            elseif itask == 2, traintrl = traintrl2; %diff task
            else               traintrl = trainall; end %both
            ntrain = length(traintrl);
            ntest  = length(testtrl);
            
            % get means from training data
            trnmuerf  = NaN(nsens,nChan,ntime,'single');
            trnNotNaN = true(nChan,1);
            mcur      = mean(dat(traintrl,:,:),1);
            sdcur     = std(dat(traintrl,:,:),[],1);
            for iang = 1:nChan
                iitrl           = traintrl(C(iang,traintrl)'==1);
                trnNotNaN(iang,1)  = nnz(iitrl)>0;
                if trnNotNaN(iang)
                    curdat          = bsxfun(@rdivide,bsxfun(@minus,dat(iitrl,:,:),mcur),sdcur);
                    trnmuerf(:,iang,:) = mean(curdat,1);
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
                    curdat          = bsxfun(@rdivide,bsxfun(@minus,dat(iitrl,:,:),mcur),sdcur);
                    muerf(:,iang,:) = mean(curdat,1);
                end
            end
            
            % project test data into task subspace from training data
            for it = 1:ntime
                coeff   = pca(trnmuerf(:,trnNotNaN,it)');                
                lo_proj = muerf(:,NotNaN,it)'*coeff(:,1:2);
                
                %get eigenvalues of covariance matrix of full data and
                %subspace projection
                [~,d]   = eig(cov(muerf(:,NotNaN,it)'));
                fullvar = sum(diag(d));
                
                [~,d]   = eig(cov(lo_proj));
                reducedvar = sum(diag(d));
                cur_variance_proj(:,it,itask,ifold,1) = [fullvar reducedvar];
                
                %get variance of full data and
                %subspace projection
                curdat  = muerf(:,NotNaN,it)';
                curmu   = mean(curdat,1);
                fullvar = sum(pdist2(curdat,curmu).^2,1)./(nChan-1);
                
                curdat     = lo_proj;
                curmu      = mean(curdat,1);
                reducedvar = sum(pdist2(curdat,curmu).^2,1)./(nChan-1);                
                cur_variance_proj(:,it,itask,ifold,2) = [fullvar reducedvar];
            end
        end
        taskvec(ifold,1) = testtask;
    end
    cfg.muvec = muvec;
    cfg.expl_var = expl_var;
    
    variance_proj = [];
    variance_proj(:,:,:,1,:) = nanmean(cur_variance_proj(:,:,:,taskvec == 0,:),4);
    variance_proj(:,:,:,2,:) = nanmean(cur_variance_proj(:,:,:,taskvec == 1,:),4);
end

if verbose
   fprintf(' - done!\n'); 
end