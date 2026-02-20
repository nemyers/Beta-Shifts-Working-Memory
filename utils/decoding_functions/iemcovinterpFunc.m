function [C,Cproj] = iemcovinterpFunc(dat,theta,thetabasis,foldvar,covtime,thetabasis2)
% [C,Cproj] = iemcovinterpFunc(dat,theta,thetabasis,foldvar,covtime,thetabasis2)    
% dat is channelXtrialXtime, theta is angles (from -pi to pi)
   
thetabasis = thetabasis(:); %make sure the basis set locations are in a column vector
nChan = length(thetabasis); %number of virtual channels of output
nChan2 = length(thetabasis2);
nsens = size(dat,2); %number of M/EEG sensors
ntrl  = size(dat,1);        
ntime = size(dat,3);
%dat   = permute(dat,[2 1 3]); %make Trial the first dimension of dat
     
if nargin < 5
    covtime = ones(1,ntime);
end
covtime = covtime(:)'; % make sure covtime is a row vector so that it works properly in the for loop below
if nargin < 4
    foldvar = [1:ntrl]';
end
if nargin < 3
    thetabasis = (-8:7)'/8 *pi;%default: 16 virtual channels from -pi to pi
end
% set up reverse encoding model parameters (i.e., design matrix)
cosfun    = @(theta,mu,pow)((0.5 + 0.5.*cos((theta-mu))).^pow);
X         = zeros(nChan+1,ntrl);
pow       = nChan - 1;
for iChan = 1:nChan
    X(iChan,:) = cosfun(theta,thetabasis(iChan),pow);
end        
X(nChan+1,:) = 1;%add constant regressor to model        

folds   = unique(foldvar);
nfold   = length(folds);

% fit model
C2    = nan(nChan+1,ntrl,ntime);
for ifold = 1:nfold
    traintrl  = find(foldvar ~= folds(ifold));
    ntraintrl = nnz(traintrl);
    testtrl   = find(foldvar == folds(ifold));
    

    %model fitting step on training data
    X1      = X(:,traintrl); %nChan+1XnTrain
    X1inv   = X1'*pinv(X1*X1');
    X1trans = X1'; %nTrainXnChan+1
    Wall    = permute(reshape(reshape(dat(traintrl,:,:),[ntraintrl nsens*ntime])'*X1inv, ...
                     [nsens ntime nChan+1]),[3 1 2]);
    clear X1inv
    
    %get residuals on training data
    res = zeros([ntraintrl nsens ntime],'single');
    for itrn = find(covtime)
       W             = Wall(:,:,itrn); %nChan+1Xnsens
       res(:,:,itrn) = X1trans*W - dat(traintrl,:,itrn); %nTrainXnsens
    end
    clear X1trans
    
    %use residuals to calculate covariance matrix
    res      = mean(res(:,:,covtime),3);
    sigma    = cov(res);
    sigmainv = pinv(sigma);
    
    %apply model to test data
    for itrn = 1:ntime
        W                  = Wall(:,:,itrn);%
        %multiply test data by weights from training data, normalized by covariance over sensors
        C2(:,testtrl,itrn) = pinv(W*W')*W*sigmainv*dat(testtrl,:,itrn)'; 
    end
end

% re-center model response (if we are only planning to use the projection
% as a summary statistic, we can think about skipping this step in the
% future)
recenter_curve = true;
C        = zeros(nChan2,ntrl,ntime,'single');
chanlist = 1:nChan2;
C2tmp    = C;
C2tmp(1:2:nChan2,:,:) = C2(1:nChan,:,:);
C2tmp(2:2:nChan2,:,:) = [(C2(1:nChan,:,:)+C2([2:nChan 1],:,:))./2];
C2 = C2tmp; clear C2tmp
if recenter_curve    
    for itrl = 1:ntrl
        iang     = theta(itrl);
        idst     = circ_dist(thetabasis2,iang);
        [~,imin] = min(abs(idst));
        shiftval = -(ceil(nChan2/2) - imin + 1);
        if shiftval<0, shiftval=shiftval+nChan2; end            
        shiftind     = [chanlist(shiftval+1:nChan2) chanlist(1:shiftval)];
        C(:,itrl,:) = squeeze(C2(shiftind,itrl,:));
    end
else
    C = C2; clear C2;
end

% project model response at each virtual sensor onto the unit vector at the correct angle 
% and sum across virtual sensors
get_projection = false;
if get_projection
    if recenter_curve
        costheta = repmat(cos(thetabasis - angspace2(ceil(nChan2/2)+1)),[1 ntrl ntime]); %get weight for each virtual channel
        Cproj    = squeeze(sum(costheta.*C,1)); %get cosine-weighted average of virtual channel responses
    else
        costheta = repmat(cos(repmat(thetabasis,[1 ntrl])-repmat(theta',[nChan2 1])),[1 1 ntime]); %get weight for each virtual channel
        Cproj    = squeeze(sum(costheta.*C,1)); %get cosine-weighted average of virtual channel responses
    end
else
    Cproj = [];
end

end