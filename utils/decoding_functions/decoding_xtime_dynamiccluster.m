%% do stats on dynamic coding
data   = decodingmatrix; %this is the group-level cross-temporal decoding matrix (subjects X time X time)
data   = (data + permute(data,[1 3 2]))./2; %make matrix symmetric
nsubs  = size(data,1);
ntimes = size(data,2);
xdiag  = nan(nsubs,ntimes,ntimes);
ydiag  = nan(nsubs,ntimes,ntimes);

for itime = 1:ntimes
    %get on-diagonal decoding value
    currdat =  data(:,itime,itime);
    %create a matrix with the on-diagonal element along all points in a row
    %(i.e. time axis 1)
    xdiag(:,itime,:) = repmat(currdat,[1 ntimes]);
    %create a matrix with the on-diagonal element along all points in a
    %column (i.e. time axis 2)
    ydiag(:,:,itime) = repmat(currdat,[1 ntimes]);
end

%subtract on-diagonal from off-diagonal values along both time axes
xdiff = data - xdiag;
ydiff = data - ydiag;

%because matrix is symmetrical, only keep the upper triangular part of it
for isub = 1:nsubs
    xdiff(isub,:,:) = triu(squeeze(xdiff(isub,:,:)));
    ydiff(isub,:,:) = triu(squeeze(ydiff(isub,:,:)));
end
%% Proper cluster testing for dynamics
[pclust,praw] = ClusterCorrectionConjunction(xdiff,ydiff,10000,0.05,0.05,false);
% print info for the largest clusters
pclustu = unique(pclust);
npclust = nnz(pclustu < 0.5);
xtlist  = repmat(tlist,[1 length(tlist)]); %tlist is the time vector, xtlist repeats tlist to create a matrix
ytlist  = repmat(tlist',[length(tlist) 1]);

for ipclust = 1:npclust
    currind  = pclust == pclustu(ipclust);
    curxrange = xtlist(permute(currind,[2 3 1]));
    curyrange = ytlist(permute(currind,[2 3 1]));
    fprintf('Cluster %02d (%03d pixels): p=%g, %d-%dms (time 1), %d-%dms (time 2).\n',ipclust,nnz(currind),pclustu(ipclust),...
            round(min(curxrange)*1000),round(max(curxrange)*1000),round(min(curyrange)*1000),round(max(curyrange)*1000));    
end