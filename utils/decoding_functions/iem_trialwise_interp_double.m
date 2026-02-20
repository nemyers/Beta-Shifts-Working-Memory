function [C,C2] = iem_trialwise_interp_double(data,angs,angs2,foldvar,angspace,angspace_targ,angspace2)
% [C,C2] = iem_trialwise_interp_double(data,angs,angs2,foldvar,angspace,angspace_targ,angspace2)
        
nChan  = length(angspace);
nChan2 = length(angspace2);
nChan_targ = length(angspace_targ);
ntime = size(data,3);
nsens = size(data,2);
ntrl  = size(data,1);

% set up reverse encoding model parameters (i.e., design matrix)
cosfun    = @(theta,mu,pow)((0.5 + 0.5.*cos((theta-mu))).^pow);
C         = zeros(nChan+nChan2+1,ntrl);
pow       = nChan - 1;
for iChan = 1:nChan
    C(iChan,:) = cosfun(angs,angspace(iChan),pow);
end        
pow       = nChan2 - 1;
for iChan = 1:nChan2
    C(iChan+nChan,:) = cosfun(angs2,angspace2(iChan),pow);
end
C(nChan+nChan2+1,:) = 1;%add constant regressor to model        

folds   = unique(foldvar);
nfold   = length(folds);

% fit model
C2    = nan(nChan+nChan2+1,ntrl,ntime);
for ifold = 1:nfold
    traintrl = find(foldvar ~= folds(ifold));
    testtrl  = find(foldvar == folds(ifold));

    C1    = C(:,traintrl);
    C1inv = C1'*pinv(C1*C1');
    Wall  = permute(reshape(reshape(data(traintrl,:,:),[nnz(traintrl) nsens*ntime])'*C1inv,[nsens ntime nChan+nChan2+1]),[3 1 2]);
    for itrn = 1:ntime
        W                  = Wall(:,:,itrn);%
        C2(:,testtrl,itrn) = pinv(W*W')*W*data(testtrl,:,itrn)'; %multiply test data by weights from training data
    end
end

% re-center model responses (with interpolation)
C         = zeros(nChan_targ,ntrl,ntime,'single');
chanlist  = 1:nChan_targ;
C3        = C2;
C2tmp     = C;
C2tmp(1:2:nChan_targ,:,:) = C2(1:nChan,:,:);
C2tmp(2:2:nChan_targ,:,:) = [(C2(1:nChan,:,:)+C2([2:nChan 1],:,:))./2];
C2 = C2tmp; clear C2tmp
for itrl = 1:ntrl
    iang     = angs(itrl);
    idst     = circ_dist(angspace_targ,iang);
    [~,imin] = min(abs(idst));
    shiftval = -(ceil(nChan_targ/2) - imin + 1);
    if shiftval<0, shiftval=shiftval+nChan_targ; end
    shiftind    = [chanlist(shiftval+1:nChan_targ) chanlist(1:shiftval)];
    C(:,itrl,:) = squeeze(C2(shiftind,itrl,:));
end

% re-center model responses
C4        = zeros(nChan2,ntrl,ntime,'single');
chanlist  = 1:nChan2;
for itrl = 1:ntrl
    iang     = angs2(itrl);
    idst     = circ_dist(angspace2,iang);
    [~,imin] = min(abs(idst));
    shiftval = -(ceil(nChan2/2) - imin + 1);
    if shiftval<0, shiftval=shiftval+nChan2; end
    shiftind     = [chanlist(shiftval+1:nChan2) chanlist(1:shiftval)];
    C4(:,itrl,:) = squeeze(C3(shiftind+nChan,itrl,:));
end
C2 = C4;
end