function [D] = mahal_trialwise(data,angs,foldvar,timelist)
% [D] = mahal_trialwise(data,angs,foldvar,timelist)

%assign stimulus variable for decoding
sortvar = angs;
trl_ind = 1:size(sortvar,1);
nbin    = 2;
pbin    = 1/4;
D       = zeros(size(data,1),nbin,size(data,3));
% time window for covariance
sigmawin = timelist >= 0 & timelist <= 0.60;

ft_progress('init','text');
th      = [];
folds   = unique(foldvar);
nfold   = length(folds);

for ifold = 1:nfold
    ft_progress(ifold/nfold);
    testtrl   = find(foldvar == folds(ifold));                
    trn_ind   = find(foldvar ~= folds(ifold));
    trn_dat   = data(trn_ind,:,:);
    X         = [cos(sortvar(trn_ind)) sin(sortvar(trn_ind))];
    X(:,end+1)= 1;
    Y         = reshape(trn_dat,[size(trn_dat,1),size(trn_dat,2)*size(trn_dat,3)]);
    b         = pinv(X)*Y;
    res       = Y-X*b;
    res       = reshape(res,[size(trn_dat,1),size(trn_dat,2),size(trn_dat,3)]);
    sigma     = covdiag(mean(res(:,:,sigmawin),3));

    %assign angle bins
    for trl = testtrl'
        tmp_angle = sortvar(trl);
        trn_angle = circ_dist(sortvar(trn_ind),tmp_angle);
        tmp_dat   = squeeze(data(trl,:,:));

        m = single([]);

        i = circ_bini(trn_angle,nbin,pbin);
        for c = 1:nbin
            m(:,c,:)   = squeeze(mean(trn_dat(i(:,c),:,:)));
            th(trl,c)  = circ_mean(trn_angle(i(:,c)));
        end

        for itrn = 1:length(timelist)                        
            if nbin < 10
                D(trl,:,itrn) = pdist2(tmp_dat(:,itrn)',m(:,:,itrn)','mahalanobis',sigma);                            
            else
                for ibin = 1:nbin
                    D(trl,ibin,itrn) = pdist([tmp_dat(:,itrn)';m(:,ibin,itrn)'],'mahalanobis',sigma);
                end
            end
        end
    end
end
cfg.th = th;
ft_progress('close');
end