function [p,praw] = ClusterCorrectionANOVA(dat, D, nSims, p_crit, p_thresh)
%% Cluster Correct across any number of factors, may only work on 1D
if nargin < 5
    p_thresh = p_crit;
end  
alpha = 0.05;
gg    = false;
[efs,F,cdfs,p]=repanova2D(dat,D,[],gg,alpha);

%[h,p,ci,stats] = ttest(dat);
pObs = p;
FObs = F;
ntest   = size(pObs,1);
maxSize = zeros(nSims,ntest);
datsize = size(dat);
%[~,permfac] = sort(rand(datsize(1),datsize(2),nSims),2);
for sim=1:nSims
    %permind = sign(rand(size(dat,1),1)-.5);
    %permind = repmat(permind,[1 datsize(2:end)]);
    simDat  = permute(shake3D(dat,2),[2 3 1]);
    [efs,F,cdfs,p]=repanova2D(simDat,D,[],gg,alpha);
    
    for itest = 1:ntest
        CC = bwconncomp(p(itest,:)<=p_thresh);
        if CC.NumObjects == 0
            maxSize(sim,itest) = 0;
        else
            tmpSize = zeros(CC.NumObjects,1);
            for c=1:CC.NumObjects
                tmpSize(c) = abs(sum(F(itest,CC.PixelIdxList{c})));
            end    
            maxSize(sim,itest) = max(tmpSize);
        end
    end
end

praw = pObs;
p    = nan(size(pObs));
for itest = 1:ntest
    CC   = bwconncomp(pObs(itest,:)<=p_crit);
    if CC.NumObjects ~= 0
        for c=1:CC.NumObjects
            ObsCluster = abs(sum(FObs(itest,CC.PixelIdxList{c})));
            p(itest,CC.PixelIdxList{c}) = (sum(maxSize(:,itest)>=ObsCluster)/nSims);
        end
    end
end
end