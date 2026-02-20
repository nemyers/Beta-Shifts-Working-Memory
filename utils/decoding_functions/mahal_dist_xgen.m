function  [mahal_dist,cfg] = mahal_dist_xgen(dat,theta,task,cfg)
% [mahal_dist,cfg] = mahal_dist_xgen(dat,theta,task,cfg)

if nargin < 3
    help mahal_dist_xgen
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
if returnConly
    mahal_dist = [];
else
    if verbose
       fprintf(' Running x-val Mahal:'); r = ''; 
    end
    folds   = unique(foldvar);
    folds   = folds(:)';
    nfold   = length(folds);
    muvec   = [];
        
    %loop over x-validation folds
    for ifold = 1:nfold
        if verbose
           m = sprintf(' Fold %d/%d',ifold,nfold); fprintf([r m]); r = repmat('\b',1,length(m)); 
        end
        testtrl    = find(foldvar == folds(ifold))';                           % get test trials
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
            
            %get training data means and covariance
            meanERP = nan(nChan,nsens,ntime,'single');
            datamu  = nan(size(dat),'single');
            for ibin = 1:nChan
                icur = traintrl(C(ibin,traintrl)'==1);
                meanERP(ibin,:,:) = mean(dat(icur,:,:),1);
                datamu(icur,:,:)  = repmat(meanERP(ibin,:,:),[nnz(icur) 1 1]);
            end
            res = dat(traintrl,:,:)-datamu(traintrl,:,:);
            sigma = covdiag(mean(res(:,:,cfg.covtime),3));
            
            %get test data means
            meanERPtst = nan(nChan,nsens,ntime,'single');
            for ibin = 1:nChan
                icur = testtrl(C(ibin,testtrl)'==1);
                meanERPtst(ibin,:,:) = mean(dat(icur,:,:),1);                
            end
            
            for itime = 1:ntime
               mahaldist = pdist2(meanERPtst(:,:,itime),meanERP(:,:,itime),'mahalanobis',sigma);
               distpair(itime,:,:,itask,ifold) = mahaldist;
            end            
        end
        taskvec(ifold,1) = testtask;
    end
    cfg.muvec  = muvec;
    mahal_dist = [];
    mahal_dist(:,:,:,:,1) = nanmean(distpair(:,:,:,:,taskvec == 0),5);
    mahal_dist(:,:,:,:,2) = nanmean(distpair(:,:,:,:,taskvec == 1),5);
end

if verbose
   fprintf(' - done!\n'); 
end