%% do stats on dynamic coding
%binmat is a [nSubjects nTimepoints nTimepoints] data matrix

%get data
data   = binmat;
data   = (data + permute(data,[1 3 2]))./2; %make data matrix symmetric
nsubs  = size(data,1);
ntimes = size(data,2);
%get diagonal values and project along either time axis
xdiag  = nan(nsubs,ntimes,ntimes);
ydiag  = nan(nsubs,ntimes,ntimes);
for itime = 1:ntimes
    currdat          =  data(:,itime,itime);
    xdiag(:,itime,:) = repmat(currdat,[1 ntimes]);
    ydiag(:,:,itime) = repmat(currdat,[1 ntimes]);
end
%get differences between diagonal and off-diagonal pattern strength
xdiff = data - xdiag; %t11 - t12
ydiff = data - ydiag; %t22 - t12

%eliminate data in lower traingular part of the matrix (because of the
%symmetry)
for isub = 1:nsubs
    xdiff(isub,:,:) = triu(squeeze(xdiff(isub,:,:)));
    ydiff(isub,:,:) = triu(squeeze(ydiff(isub,:,:)));
end
%% Proper cluster testing for dynamics
%check if off-diagonal points are significantly lower than both
%intersections with the diagonal
[pclust,praw] = ClusterCorrectionConjunction(xdiff,ydiff,10000,0.05,0.05,false);
%
plot_result   = true;
clc
pclustu = unique(pclust);
npclust = nnz(pclustu < 0.5);
%tlist is the [1 nTimepoints] time vector
xtlist  = repmat(tlist,[1 length(tlist)]);
ytlist  = repmat(tlist',[length(tlist) 1]);

if plot_result
    pmap  = zeros(size(pclust),'single');
end
for ipclust = 1:npclust
    currind  = pclust == pclustu(ipclust);
    curxrange = xtlist(permute(currind,[2 3 1]));
    curyrange = ytlist(permute(currind,[2 3 1]));
    fprintf('Cluster %02d (%03d pixels): p=%g, %d-%dms (train), %d-%dms (test).\n',ipclust,nnz(currind),pclustu(ipclust),round(min(curxrange)*1000),round(max(curxrange)*1000),round(min(curyrange)*1000),round(max(curyrange)*1000));    
    if plot_result
        pmap(currind) = ipclust;
    end
end
if plot_result & npclust > 0
    close all
    figure
    set(gcf,'color','white')    
    imagesc(tlist,tlist,permute(pmap,[2 3 1]),[-npclust npclust]);
    hold on
    plot(tlist([1 end]),[0 0],'k-')
    plot([0 0],tlist([1 end]),'k-')
    plot(tlist([1 end]),tlist([1 end]),'k-')
    axis xy
    axis square
end
%% plot cross-generalization time
close all
clc
figure
set(gcf,'color','white')

xlims = tlist([1 end]);
ylims = [-0.1 0.4];
wins  = [0.05 0.075];
nwins = size(wins,1);
cols  = colormap(jet(nwins));
for iwin = 1:nwins
    win = wins(iwin,:);
    iit = tlist >= win(1) & tlist <= win(2);
    dat = squeeze(mean(mean(data(:,iit,:),2),1));
    xp  = [win win(2:-1:1)];
    yp  = sort([ylims ylims]);
    cp  = cols(iwin,:);
    patch(xp,yp,cp,'Edgecolor','none','FaceAlpha',0.25);
    hold on
    plot(tlist,dat,'k-','color',cp,'linewidth',2);    
end