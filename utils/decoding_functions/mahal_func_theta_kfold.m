function  [distance_cos,distance_difference,pred_cond,distances] = ...
        mahal_func_theta_kfold(data,conditions,theta,n_folds)

%%
train_partitions = cvpartition(conditions,'KFold',n_folds); % split data n times using Kfold
distances=nan(length(unique(conditions)),size(data,1),size(data,3)); % prepare for output 
distance_difference=nan(size(data,1),size(data,3));

theta=circ_dist(theta,0);

u_theta=unique(theta);

theta_dist=circ_dist2(u_theta',theta)';

uc = unique(conditions);
nc = length(uc);

for tst=1:n_folds % run for each fold
    trn_ind = training(train_partitions,tst); % get training trial rows
    tst_ind = test(train_partitions,tst); % get test trial rows
    trn_dat = data(trn_ind,:,:); % isolate training data
    tst_dat = data(tst_ind,:,:); % isolate test data
    trn_cond =conditions(trn_ind);
    tst_cond =conditions(tst_ind);
    
    m=double(nan(length(unique(conditions)),size(data,2),size(data,3)));

    % average trials over each condition of training data
    for c=1:nc
        m(c,:,:)=mean(trn_dat(trn_cond==uc(c),:,:),1);  
    end
    for t=1:size(data,3) % decode at each time-point
        if ~isnan(trn_dat(:,:,t))
            % compute pair-wise mahalabonis distance between test-trials
            % and averaged training data, using the covariance matrix
            % computed from the training data
            temp=pdist2(squeeze(m(:,:,t)), squeeze(tst_dat(:,:,t)),'mahalanobis',covdiag(trn_dat(:,:,t)));            
            distances(:,tst_ind,t)=temp;
        else
            distances(:,tst_ind,t)=nan;
        end
    end
end

distance_cos=-mean(bsxfun(@times,cos(theta_dist)',distances),1);

if nargout > 3
    [~,I]=(min(distances,[],1));
    pred_cond=squeeze(I);
end

if nargout > 1
    for c=1:length(unique(conditions))
        temp=round(circ_dist(u_theta(c),u_theta),4);
        temp(temp==round(pi,4))=round(-pi,4);
        [~,i]=sort(temp);
        distances(1:length(unique(conditions)),conditions==c,:)=distances(i,conditions==c,:);

        ind=find(conditions==c);

        distance_difference(ind,:)=mean(distances(setdiff(unique(conditions),conditions(c)),...
            ind,:),1)-distances(c,ind,:);
    end
end

if nargout > 2
    distances=-bsxfun(@minus,distances,mean(distances,1));
end
end