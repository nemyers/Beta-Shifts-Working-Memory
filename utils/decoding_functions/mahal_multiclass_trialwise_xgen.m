function  [C2,cfg] = mahal_multiclass_trialwise_xgen(dat,labels,task,cfg)
% [C2,cfg] = mahal_multiclass_trialwise(dat,labels,task,cfg)

% authors: M Stokes, E Spaak, N Myers

if nargin < 2
    help mahal_multiclass_trialwise_xgen
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

%flag setting whether covariance is estimated only on same-task trials (and
%applied to other-task trials also)
if ~isfield(cfg,'covwithintask') || isempty(cfg.covwithintask)
    cfg.covwithintask = false; end
covwithintask    = cfg.covwithintask;

%flag setting whether covariance is estimated over all training trials (and
%applied to both same and other-task trials)
if ~isfield(cfg,'covavg') || isempty(cfg.covavg)
    cfg.covavg = false; end
covavg = cfg.covavg;

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

%return only design matrix and quit?
if ~isfield(cfg,'returnConly') || isempty(cfg.returnConly)
    cfg.returnConly = false; end
returnConly = cfg.returnConly;
    
%return only covariance matrix for each fold and quit?
if ~isfield(cfg,'getSigmaOnly') || isempty(cfg.getSigmaOnly)
    cfg.getSigmaOnly = false; end
getSigmaOnly = cfg.getSigmaOnly;

%match distributions of training trials in terms of distance to test
%trials?
if ~isfield(cfg,'match_train_trial_distance_to_test_trial') || isempty(cfg.match_train_trial_distance_to_test_trial)
    cfg.match_train_trial_distance_to_test_trial = true; end
match_train_trial_distance_to_test_trial = cfg.match_train_trial_distance_to_test_trial;

%match numbers of training trials?
if ~isfield(cfg,'match_training_trial_numbers') || isempty(cfg.match_training_trial_numbers)
    cfg.match_training_trial_numbers = true; end
match_training_trial_numbers = cfg.match_training_trial_numbers;

%demean data? (instead of adding a constant regressor)
if ~isfield(cfg,'demean') || isempty(cfg.demean)
    cfg.demean = false; end
demean = cfg.demean;
if demean, cfg.add_constant = false; end

%check if cross-task decoding needs to be flipped?
if ~isfield(cfg,'xvalflip') || isempty(cfg.xvalflip)
    cfg.xvalflip = false; end
xvalflip = cfg.xvalflip;

%add time window to determine cross-task decoding flip?
if ~isfield(cfg,'xvalflipwin') || isempty(cfg.xvalflipwin)
    cfg.xvalflipwin = true(ntime,1); end
xvalflipwin = cfg.xvalflipwin;

%add time window to determine cross-task decoding flip?
if ~isfield(cfg,'nuisance') || isempty(cfg.nuisance)
    cfg.nuisance = [];
    cfg.add_nuisance = false; end
nuisance = cfg.nuisance;
add_nuisance = cfg.add_nuisance;

X = dummyvar(labels);
X(:,all(X==repmat(mean(X,1),[size(X,1) 1]),1)) = [];
nstim = size(X,2);
X = [X ones(length(labels),add_constant)];
if add_nuisance
    nnuis = size(nuisance,2);
    X = [X nuisance]; end
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
            
            testtrl_all  = find(foldvar == folds(ifold))';% get test trials
            testtasks    = unique(task(testtrl_all));
            ntesttask    = length(testtasks);
        
            for itesttask = 1:ntesttask
                testtrl = testtrl_all(task(testtrl_all)==testtasks(itesttask));
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
                    %diffthresh   = 5;
                    diffthresh   = 0;
                    iter         = 1;
                    maxreduction = 0.20; %maximum proportion of trials to throw out
                    %while abs(diff(mu))>diffthresh
                    while diff(mu)>diffthresh
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

                testtrl = testtrl';

                %loop over tasks
                for itask = 1:2
                    if itask == 1, traintrl = traintrl1; else traintrl = traintrl2; end
                    %get condition means for training data
                    if add_nuisance
                        Xtmp = X(traintrl,nstim+add_constant+1:end);
                        Xtmpmu = mean(Xtmp,1);
                        Xtmp = bsxfun(@minus,Xtmp,Xtmpmu);
                        [coeff score] = pca(Xtmp);
                    end
                    
                    if add_nuisance
                        pC       = pinv([X(traintrl,1:nstim+add_constant) score(:,1)]);
                    else
                        pC       = pinv([X(traintrl,1:nstim+add_constant)]);
                    end
                    dtmp     = dat(traintrl,:,:);
                    if demean, Dmu = mean(dtmp,1); dtmp = bsxfun(@minus,dtmp,Dmu); end
                    dtmp     = reshape(dtmp,[nnz(traintrl) nsens*ntime]);
                    b        = pC*dtmp;
                    %res      = dat(traintrl,:,:) - reshape([b'*X(traintrl,:)']',[nnz(traintrl) nsens ntime]);
                    D        = reshape(b(1:nstim+add_constant,:),[nstim+add_constant nsens ntime]);
                    if add_constant, Dmu = D(end,:,:); D(end,:,:) = []; end

                    % get test data
                    dattest = dat(testtrl,:,:);
                    if add_constant || demean
                        dattest = bsxfun(@minus,dattest,Dmu); end
                    if add_nuisance
                        Xtesttmp = X(testtrl,nstim+add_constant+1:end);
                        Xtesttmp = bsxfun(@minus,Xtesttmp,Xtmpmu);
                        Xtesttmp = Xtesttmp*coeff(:,1);
                        dattest = dattest - reshape(Xtesttmp*b(nstim+add_constant+1:end,:),[length(testtrl) nsens ntime]); end

                    %sigma    = covdiag(mean(res(:,:,covtime),3)); %covdiag or cov
                    if add_nuisance
                        if (~covbytime | itask==1 | (covbytime & ~covwithintask)), res = dat(traintrl,:,:) - reshape([b'*[X(traintrl,1:nstim+add_constant) score(:,1)]']',[nnz(traintrl) nsens ntime]); end;
                    else
                        if (~covbytime | itask==1 | (covbytime & ~covwithintask)), res = dat(traintrl,:,:) - reshape([b'*X(traintrl,1:nstim+add_constant)']',[nnz(traintrl) nsens ntime]); end;
                    end
                    if (~covbytime & (itask==1 | ~covwithintask)), sigma = covdiag(mean(res(:,:,covtime),3)); end %covdiag or cov                    
                    
                    for ti = 1:ntime
                        if covbytime
                            sigma = covdiag(res(:,:,ti)); end                            
                        C2(:,testtrl,ti,itask) = pdist2(dattest(:,:,ti),D(:,:,ti),'mahalanobis',sigma)';
                    end
                    
                    if xvalflip & itask == 2
                        dattest = dat(traintrl1,:,:);
                        if add_constant || demean
                            dattest = bsxfun(@minus,dattest,Dmu); end
                        if add_nuisance
                            Xtesttmp = X(traintrl1,nstim+add_constant+1:end);
                            Xtesttmp = bsxfun(@minus,Xtesttmp,Xtmpmu);
                            Xtesttmp = Xtesttmp*coeff(:,1);
                            dattest = dattest - reshape(Xtesttmp*b(nstim+add_constant+1:end,:),[length(traintrl1) nsens ntime]); end
                        %dattest = dattest - reshape(X(traintrl1,nstim+add_constant+1:end)*b(nstim+add_constant+1:end,:),[length(traintrl1) nsens ntime]); end
                        Ctmp = [];
                        for ti = 1:ntime
                            if covbytime
                                sigma = covdiag(res(:,:,ti)); end                            
                            Ctmp(:,:,ti) = pdist2(dattest(:,:,ti),D(:,:,ti),'mahalanobis',sigma)';
                        end
                        
                        %Cmu = [];
                        %Cmu(1,:) = mean(mean(Ctmp(:,X(traintrl1,1)==1,xvalflipwin),2),3);
                        %Cmu(2,:) = mean(mean(Ctmp(:,X(traintrl1,1)~=1,xvalflipwin),2),3);
                        %Cmu(2,:) = Cmu(2,[2 1]);
                        %Cmu = mean(Cmu,1);
                        %flipsign = sign(-diff(Cmu));
                        %cfg.flipvec(testtrl,1) = flipsign;
                        
                        Ctmp(:,X(traintrl1,1)~=1,:) = Ctmp([2 1],X(traintrl1,1)~=1,:);
                        Cmu = mean(mean(Ctmp(:,:,xvalflipwin),2),3);
                        flipsign = sign(diff(Cmu));
                        cfg.flipvec(testtrl,1) = flipsign;
                    end
                end
            end
        end
    case {'euclid' 'euclidean' 'e' 'euc'}
        for ifold = 1:nfold
            if verbose
               m = sprintf(' Fold %d/%d',ifold,nfold); fprintf([r m]); r = repmat('\b',1,length(m)); 
            end
            
            testtrl_all  = find(foldvar == folds(ifold))';% get test trials
            testtasks    = unique(task(testtrl_all));
            ntesttask    = length(testtasks);
        
            for itesttask = 1:ntesttask
                testtrl = testtrl_all(task(testtrl_all)==testtasks(itesttask));
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
                    %diffthresh   = 5;
                    diffthresh   = 0;
                    iter         = 1;
                    maxreduction = 0.20; %maximum proportion of trials to throw out
                    %while abs(diff(mu))>diffthresh
                    while diff(mu)>diffthresh
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

                testtrl = testtrl';

                %loop over tasks
                for itask = 1:2
                    if itask == 1, traintrl = traintrl1; else traintrl = traintrl2; end
            
                    %get condition means for training data
                    pC       = pinv(X(traintrl,:));
                    dtmp     = dat(traintrl,:,:);
                    if demean, Dmu = mean(dtmp,1); dtmp = bsxfun(@minus,dtmp,Dmu); end
                    dtmp     = reshape(dtmp,[nnz(traintrl) nsens*ntime]);
                    b        = pC*dtmp;
                    D        = reshape(b(1:nstim+add_constant,:),[nstim+add_constant nsens ntime]);
                    if add_constant, Dmu = D(end,:,:); D(end,:,:) = []; end
                    
                    % get test data
                    dattest = dat(testtrl,:,:);
                    if add_constant || demean
                        dattest = bsxfun(@minus,dattest,Dmu); end
                    if add_nuisance
                        dattest = dattest - reshape(X(testtrl,nstim+add_constant+1:end)*b(nstim+add_constant+1:end,:),[length(testtrl) nsens ntime]); end
                    
                    for ti = 1:ntime                   
                        C2(:,testtrl,ti,itask) = pdist2(dattest(:,:,ti),D(:,:,ti),'euclidean')';
                    end
                    
                    if xvalflip & itask == 2
                        dattest = dat(traintrl1,:,:);
                        if add_constant || demean
                            dattest = bsxfun(@minus,dattest,Dmu); end
                        if add_nuisance
                        dattest = dattest - reshape(X(testtrl,nstim+add_constant+1:end)*b(nstim+add_constant+1:end,:),[length(testtrl) nsens ntime]); end
                    
                        Ctmp = [];
                        for ti = 1:ntime      
                            Ctmp(:,:,ti) = pdist2(dattest(:,:,ti),D(:,:,ti),'euclidean')';
                        end
                        
                        Ctmp(:,X(traintrl1,1)~=1,:) = Ctmp([2 1],X(traintrl1,1)~=1,:);
                        Cmu = mean(mean(Ctmp(:,:,xvalflipwin),2),3);
                        flipsign = sign(diff(Cmu));
                        cfg.flipvec(testtrl,1) = flipsign;
                    end
                end
            end
    end
end

if verbose
   fprintf(' - done!\n'); 
end