function [C] = iem_trialwise_pca(data,angs,angspace,foldvar,npcs)
% [C] = iem_trialwise_pca(data,angs,angspace,foldvar,npcs)     
        
nChan = length(angspace);
ntime = size(data,3);
nsens = size(data,2);
ntrl  = size(data,1);        
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
    traintrl = find(foldvar ~= folds(ifold));
    testtrl  = find(foldvar == folds(ifold));
    
    zmu   = mean(data(traintrl,:,:),   1);
    zsd   = std( data(traintrl,:,:),[],1);
    zdat  = bsxfun(@rdivide,bsxfun(@minus,data,zmu),zsd);
    
    pcdat   = zeros(ntrl,npcs,ntime,'single');
    for it = 1:ntime
       [u,~,~]         = svd(cov(zdat(traintrl,:,it)));
       pcdat(:,:,it)   = zdat(:,:,it)*u(:,1:npcs);
    end
    clear zdat
    
    C1    = C(:,traintrl);
    C1inv = C1'*pinv(C1*C1');
    Wall  = permute(reshape(reshape(pcdat(traintrl,:,:),[nnz(traintrl) npcs*ntime])'*C1inv,[npcs ntime nChan+1]),[3 1 2]);
    for itrn = 1:ntime
        W                  = Wall(:,:,itrn);%
        C2(:,testtrl,itrn) = pinv(W*W')*W*pcdat(testtrl,:,itrn)'; %multiply test data by weights from training data
    end
end

% re-center model response
C        = zeros(nChan,ntrl,ntime,'single');
chanlist = 1:nChan; 
for itrl = 1:ntrl
    iang     = angs(itrl);
    idst     = circ_dist(angspace,iang);
    [~,imin] = min(abs(idst));
    shiftval = -(ceil(nChan/2) - imin + 1);
    if shiftval<0, shiftval=shiftval+nChan; end            
    shiftind     = [chanlist(shiftval+1:nChan) chanlist(1:shiftval)];
    C(:,itrl,:) = squeeze(C2(shiftind,itrl,:));
end

end