function  cos_amp = mahalTune_func_fast(data,theta,num_features,weights,foldvar)

% computes the mahalanobis distances between the data of single test-trial of a particular orientation and the rest of the data,
% averaged into specific orientation bins relative to orientation of the test-trial.

%% input
% data is trial by channel by time

% theta is angles in radians (-pi to pi)

% num_comps is the number of pca components that want to be used, default
% is all

% weights is the weight give to each trial, default is 1

%% output
% cos_amp is the cosine amplitude of the tuning curve, interpreted as
% decoding accuracy
if nargin<4||isempty(weights)
    weights = ones(size(data,1),1);
end
if nargin<3||isempty(num_features)
   num_features=size(data,2);
end
if num_features>size(data,2)
    error('number of features cannot be higher than the input');
end
%%
cos_amp=nan(size(data,1),size(data,3));
theta_mat=(circ_dist2(theta,theta)); % compute circular distances of the presented sitmuli between all trials
theta_cos_mat=cos(theta_mat); % get the cosine
theta_cos_mat(logical(eye(size(theta_cos_mat)))) = 0; % set the diagonal to zero, so that there is no bias later on

ntrl    = size(data,1);
foldmat = false(ntrl);
for itrl = 1:ntrl
   foldmat(itrl,foldvar==foldvar(itrl)) = true;
end
theta_cos_mat(foldmat) = 0;

weights = weights*weights';
for t=1:size(data,3)
    if ~isnan(data(:,:,t))
        [~, data_pca_temp] = pca(data(:,:,t));
        temp = squareform(pdist(zscore(data_pca_temp(:,1:num_features)),'euclid'));
        cos_amp(:,t)=(mean(theta_cos_mat.*-(temp).*weights,2))... %convolve the trialwise distances with the the cosine trial-wise distances
            -mean(theta_cos_mat.*repmat(mean(-temp.*weights,2),[1,size(temp,1)]),2); % subtract the mean vector means from each trial to normalize
        
        %cos_amp(:,t)=(mean(theta_cos_mat.*-temp,2));% - mean(theta_cos_mat.*repmat(mean(-temp,2),[1,size(temp,1)]),2); % subtract the mean vector means from each trial to normalize
    end
end
end


