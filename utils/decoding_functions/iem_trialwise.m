function [C thetaout] = iem_trialwise(data,angs,angspace,foldvar,exclusion,recenter)
%        [C thetaout] = iem_trialwise(data,angs,angspace,foldvar,exclusion,recenter)     
        
nChan = length(angspace);
ntrl  = size(data,1);  
nsens = size(data,2);
ntime = size(data,3);
if nargin < 5
   exclusion = zeros(ntrl,1);
end
if nargin < 6
    recenter = true;
end
% set up reverse encoding model parameters (i.e., design matrix)
cosfun    = @(theta,mu,pow)((0.5 + 0.5.*cos((theta-mu))).^pow);
C         = zeros(nChan+1,ntrl);
for iChan = 1:nChan
    pow        = nChan - 1;
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

% re-center model response
C        = zeros(nChan,ntrl,ntime,'single');
chanlist = 1:nChan; 
if recenter
    thetaout = angs;
    for itrl = 1:ntrl
        iang     = angs(itrl);
        idst     = circ_dist(angspace,iang);
        [~,imin] = min(abs(idst));
        thetaout(itrl) = idst(imin);
        shiftval = -(ceil(nChan/2) - imin + 1);
        if shiftval<0, shiftval=shiftval+nChan; end            
        shiftind     = [chanlist(shiftval+1:nChan) chanlist(1:shiftval)];
        C(:,itrl,:) = squeeze(C2(shiftind,itrl,:));
    end
else
    thetaout = angs;
    C = C2(chanlist,:,:);
end

end