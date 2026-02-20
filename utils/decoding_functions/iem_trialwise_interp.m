function [C] = iem_trialwise_interp(data,angs,angspace,foldvar,angspace2,exclusion,recenter)
% [C] = iem_trialwise_interp(data,angs,angspace,foldvar,angspace2,exclusion,recenter)     
        
nChan = length(angspace);
nChan2 = length(angspace2);
ntime = size(data,3);
nsens = size(data,2);
ntrl  = size(data,1);        

if nargin < 6
   exclusion = zeros(ntrl,1);
end
if nargin < 7
    recenter = true;
end
% set up reverse encoding model parameters (i.e., design matrix)
cosfun    = @(theta,mu,pow)((0.5 + 0.5.*cos((theta-mu))).^pow);
C         = zeros(nChan+1,ntrl);
pow       = nChan - 1;
for iChan = 1:nChan
    C(iChan,:) = cosfun(angs,angspace(iChan),pow);
end        
C(nChan+1,:) = 1;%add constant regressor to model        

folds   = unique(foldvar);
nfold   = length(folds);

% fit model
C2    = nan(nChan+1,ntrl,ntime);
for ifold = 1:nfold
    traintrl = find(foldvar ~= folds(ifold) & ~exclusion);
    testtrl  = find(foldvar == folds(ifold));

    C1    = C(:,traintrl);
    C1inv = C1'*pinv(C1*C1');
    Wall  = permute(reshape(reshape(data(traintrl,:,:),[nnz(traintrl) nsens*ntime])'*C1inv,[nsens ntime nChan+1]),[3 1 2]);
    for itrn = 1:ntime
        W                  = Wall(:,:,itrn);%
        C2(:,testtrl,itrn) = pinv(W*W')*W*data(testtrl,:,itrn)'; %multiply test data by weights from training data
    end
end

% re-center model responses (with interpolation)
C         = zeros(nChan2,ntrl,ntime,'single');
chanlist  = 1:nChan2;
C2tmp     = C;
C2tmp(1:2:nChan2,:,:) = C2(1:nChan,:,:);
C2tmp(2:2:nChan2,:,:) = [(C2(1:nChan,:,:)+C2([2:nChan 1],:,:))./2];
C2 = C2tmp; clear C2tmp
if recenter
    for itrl = 1:ntrl
        iang     = angs(itrl);
        idst     = circ_dist(angspace2,iang);
        [~,imin] = min(abs(idst));
        shiftval = -(ceil(nChan2/2) - imin + 1);
        if shiftval<0, shiftval=shiftval+nChan2; end
        shiftind    = [chanlist(shiftval+1:nChan2) chanlist(1:shiftval)];
        C(:,itrl,:) = squeeze(C2(shiftind,itrl,:));
    end
else
   C = C2(chanlist,:,:); 
end
end