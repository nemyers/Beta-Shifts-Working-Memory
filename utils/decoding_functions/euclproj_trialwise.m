function [D] = euclproj_trialwise(data,angs,foldvar,timelist)
% [D] = euclproj_trialwise(data,angs,foldvar,timelist)

%assign stimulus variable for decoding
sortvar = angs;
nbin    = 2;
pbin    = 1/4;
D       = zeros(size(data,1),size(data,3));

ft_progress('init','text');
th      = [];
folds   = unique(foldvar);
nfold   = length(folds);

for ifold = 1:nfold
    ft_progress(ifold/nfold);
    testtrl   = find(foldvar == folds(ifold));                
    trn_ind   = find(foldvar ~= folds(ifold));
    trn_dat   = data(trn_ind,:,:);
    
    %assign angle bins
    for trl = testtrl'
        tmp_angle = sortvar(trl);
        trn_angle = circ_dist(sortvar(trn_ind),tmp_angle);
        tmp_dat   = squeeze(data(trl,:,:));

        m = single([]);
        i = circ_bini(trn_angle,nbin,pbin);
        for c = 1:nbin
            m(:,c,:)   = squeeze(mean(trn_dat(i(:,c),:,:)));
        end
        ms = squeeze(diff(m,[],2));
        
        for itrn = 1:length(timelist)
            D(trl,itrn) = ms(:,itrn)'*tmp_dat(:,itrn);
        end
    end
end
ft_progress('close');
end